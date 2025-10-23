"""
Komega CG Real Module - Real CG solver implementation

This module corresponds to the Fortran komega_cg_r module and provides
the CG (Conjugate Gradient) solver for real linear systems.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

from typing import List, Optional, Tuple

import numpy as np

from komega_math import dabsmax, daxpy, dcopy, ddotMPI, dscal, get_global_math
from komega_parameter import get_global_params, initialize_global_params
from komega_vals_r import cleanup_vals_r, get_global_vals_r, initialize_vals_r
from komega_vecs_r import cleanup_vecs_r, get_global_vecs_r, initialize_vecs_r


class KomegaCGR:
    """
    CG (Conjugate Gradient) solver for real linear systems.

    This class implements the CG algorithm with shifted equations
    for solving multiple frequency points simultaneously.
    """

    def __init__(self):
        """Initialize the CG solver."""
        self.params = get_global_params()
        self.math = get_global_math()
        self.vals_r = get_global_vals_r()
        self.vecs_r = get_global_vecs_r()

        # Status tracking
        self.initialized = False

    def init(
        self,
        ndim: int,
        nl: int,
        nz: int,
        z: np.ndarray,
        itermax: int,
        threshold: float,
        comm: Optional[int] = None,
    ) -> np.ndarray:
        """
        Initialize the CG solver.

        Parameters
        ----------
        ndim : int
            Dimension of the Hamiltonian matrix
        nl : int
            Dimension of the projection space
        nz : int
            Number of frequency points
        z : np.ndarray
            Frequency array (real)
        itermax : int
            Maximum number of iterations
        threshold : float
            Convergence threshold
        comm : int, optional
            MPI communicator (for MPI version)

        Returns
        -------
        np.ndarray
            Initialized solution array (nl x nz)
        """
        # Initialize global parameters
        initialize_global_params(ndim, nl, nz, itermax, threshold, comm)

        # Initialize storage modules
        initialize_vals_r(z, itermax)
        initialize_vecs_r(ndim, nl, nz, itermax)

        # Initialize solution array
        x = np.zeros((nl, nz), dtype=float)

        self.initialized = True
        return x

    def shifted_equation(self, r_l: np.ndarray, x: np.ndarray) -> None:
        """
        Update shifted equations for all frequencies.

        Parameters
        ----------
        r_l : np.ndarray
            Local residual vector
        x : np.ndarray
            Solution array (modified in place)
        """
        if not self.initialized:
            return

        for iz in range(self.params.nz):
            if self.params.lz_conv[iz]:
                continue

            # Update π values
            pi_new = (
                1.0 + self.vals_r.alpha * (self.vals_r.z[iz] - self.vals_r.z_seed)
            ) * self.vals_r.pi[
                iz
            ] - self.vals_r.alpha * self.vals_r.beta / self.vals_r.alpha_old * (
                self.vals_r.pi_old[iz] - self.vals_r.pi[iz]
            )

            # Update search direction
            self.vecs_r.p[:, iz] = (
                r_l / self.vals_r.pi[iz]
                + (self.vals_r.pi_old[iz] / self.vals_r.pi[iz]) ** 2
                * self.vals_r.beta
                * self.vecs_r.p[:, iz]
            )

            # Update solution
            daxpy(
                self.vals_r.pi[iz] / pi_new * self.vals_r.alpha,
                self.vecs_r.p[:, iz],
                x[:, iz],
            )

            # Update π values
            self.vals_r.pi_old[iz] = self.vals_r.pi[iz]
            self.vals_r.pi[iz] = pi_new

            # Save for restart if needed
            if self.params.itermax > 0:
                self.vals_r.pi_save[iz, self.params.iter] = pi_new

    def seed_switch(self, v2: np.ndarray, status: List[int]) -> None:
        """
        Perform seed switching to find the best frequency.

        Parameters
        ----------
        v2 : np.ndarray
            Residual vector
        status : List[int]
            Status array (modified in place)
        """
        if not self.initialized:
            return

        # Find minimum |π| among unconverged frequencies
        unconverged_indices = [
            i for i, conv in enumerate(self.params.lz_conv) if not conv
        ]
        if not unconverged_indices:
            status[2] = 0
            return

        pi_abs = np.abs(self.vals_r.pi[unconverged_indices])
        min_idx = np.argmin(pi_abs)
        status[2] = unconverged_indices[min_idx]

        if abs(self.vals_r.pi[status[2]]) < self.params.almost0:
            status[1] = 3

        if status[2] != self.params.iz_seed:
            # Update seed frequency
            self.params.iz_seed = status[2]
            self.vals_r.z_seed = self.vals_r.z[self.params.iz_seed]

            # Update α and ρ
            self.vals_r.alpha *= (
                self.vals_r.pi_old[self.params.iz_seed]
                / self.vals_r.pi[self.params.iz_seed]
            )
            self.vals_r.rho /= self.vals_r.pi_old[self.params.iz_seed] ** 2

            # Scale vectors
            scale = 1.0 / self.vals_r.pi[self.params.iz_seed]
            dscal(scale, v2)
            self.vals_r.scale_pi_values(scale)

            # Scale old vectors
            scale_old = 1.0 / self.vals_r.pi_old[self.params.iz_seed]
            dscal(scale_old, self.vecs_r.v3)
            self.vals_r.pi_old *= scale_old

    def update(
        self,
        v12: np.ndarray,
        v2: np.ndarray,
        x: np.ndarray,
        r_l: np.ndarray,
        status: List[int],
    ) -> None:
        """
        Update the CG iteration.

        Parameters
        ----------
        v12 : np.ndarray
            Working vector 12
        v2 : np.ndarray
            Residual vector
        x : np.ndarray
            Solution array (modified in place)
        r_l : np.ndarray
            Local residual vector
        status : List[int]
            Status array (modified in place)
        """
        if not self.initialized:
            return

        # Increment iteration counter
        self.params.increment_iteration()
        status[:] = [0, 0, 0]

        # Update ρ and β
        rho_old = self.vals_r.rho
        self.vals_r.rho = ddotMPI(v2, v2)

        if self.params.iter == 1:
            self.vals_r.beta = 0.0
        else:
            self.vals_r.beta = self.vals_r.rho / rho_old

        # Update working vector
        v12[:] = self.vals_r.z_seed * v2 - v12

        # Update α
        self.vals_r.alpha_old = self.vals_r.alpha
        alpha_denom = (
            ddotMPI(v2, v12) - self.vals_r.beta * self.vals_r.rho / self.vals_r.alpha
        )

        if abs(alpha_denom) < self.params.almost0:
            status[1] = 2

        self.vals_r.alpha = self.vals_r.rho / alpha_denom

        # Save for restart
        if self.params.itermax > 0:
            self.vals_r.alpha_save[self.params.iter - 1] = self.vals_r.alpha
            self.vals_r.beta_save[self.params.iter - 1] = self.vals_r.beta
            self.vecs_r.save_r_l(r_l, self.params.iter)

        # Update shifted equations
        self.shifted_equation(r_l, x)

        # Update residual vector
        v12[:] = (
            (1.0 + self.vals_r.alpha * self.vals_r.beta / self.vals_r.alpha_old) * v2
            - self.vals_r.alpha * v12
            - self.vals_r.alpha
            * self.vals_r.beta
            / self.vals_r.alpha_old
            * self.vecs_r.v3
        )
        dcopy(v2, self.vecs_r.v3)
        dcopy(v12, v2)

        # Perform seed switching
        self.seed_switch(v2, status)

        # Check convergence
        v12[0] = np.sqrt(ddotMPI(v2, v2))
        self.params.resnorm = v12[0]

        # Update convergence flags
        for iz in range(self.params.nz):
            if abs(v12[0] / self.vals_r.pi[iz]) < self.params.threshold:
                self.params.lz_conv[iz] = True

        # Set status
        if self.params.resnorm < self.params.threshold:
            status[0] = -self.params.iter
            status[1] = 0
        elif self.params.iter == self.params.itermax:
            status[0] = -self.params.iter
            status[1] = 1
        elif status[1] == 2:
            status[0] = -self.params.iter
        elif status[1] == 3:
            status[0] = -self.params.iter
        else:
            status[0] = self.params.iter
            status[1] = 0

    def get_coefficients(self) -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
        """
        Get saved coefficients for restart.

        Returns
        -------
        tuple
            (alpha_save, beta_save, z_seed, r_l_save)
        """
        if not self.initialized:
            return np.array([]), np.array([]), 0.0, np.array([])

        alpha_save, beta_save, _ = self.vals_r.get_saved_values(self.params.iter)
        r_l_save = self.vecs_r.get_saved_r_l(self.params.iter)

        return alpha_save, beta_save, self.vals_r.z_seed, r_l_save

    def get_vectors(self) -> np.ndarray:
        """
        Get old residual vector.

        Returns
        -------
        np.ndarray
            Old residual vector
        """
        if not self.initialized:
            return np.array([])

        return self.vecs_r.get_v3()

    def get_residual(self) -> np.ndarray:
        """
        Get residual norms for all frequencies.

        Returns
        -------
        np.ndarray
            Residual norms
        """
        if not self.initialized:
            return np.array([])

        return self.params.resnorm / np.abs(self.vals_r.pi)

    def finalize(self) -> None:
        """Finalize the CG solver and clean up memory."""
        if self.initialized:
            cleanup_vals_r()
            cleanup_vecs_r()
            self.initialized = False


# Global CG solver instance
_global_cg_r = KomegaCGR()


def get_global_cg_r() -> KomegaCGR:
    """
    Get the global CG solver instance.

    Returns
    -------
    KomegaCGR
        Global CG solver instance
    """
    return _global_cg_r


# Convenience functions for direct access
def komega_CG_R_init(
    ndim: int,
    nl: int,
    nz: int,
    z: np.ndarray,
    itermax: int,
    threshold: float,
    comm: Optional[int] = None,
) -> np.ndarray:
    """Initialize the CG solver."""
    return _global_cg_r.init(ndim, nl, nz, z, itermax, threshold, comm)


def komega_CG_R_update(
    v12: np.ndarray, v2: np.ndarray, x: np.ndarray, r_l: np.ndarray, status: List[int]
) -> None:
    """Update the CG iteration."""
    _global_cg_r.update(v12, v2, x, r_l, status)


def komega_CG_R_getcoef() -> Tuple[np.ndarray, np.ndarray, float, np.ndarray]:
    """Get saved coefficients for restart."""
    return _global_cg_r.get_coefficients()


def komega_CG_R_getvec() -> np.ndarray:
    """Get old residual vector."""
    return _global_cg_r.get_vectors()


def komega_CG_R_getresidual() -> np.ndarray:
    """Get residual norms for all frequencies."""
    return _global_cg_r.get_residual()


def komega_CG_R_finalize() -> None:
    """Finalize the CG solver."""
    _global_cg_r.finalize()
