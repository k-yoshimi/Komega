"""
Komega COCG Module - COCG solver implementation

This module corresponds to the Fortran komega_cocg module and provides
the COCG (Conjugate Orthogonal Conjugate Gradient) solver for complex linear systems.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

import numpy as np
from typing import Tuple, Optional, List
from .komega_parameter import get_global_params, initialize_global_params
from .komega_math import get_global_math, zdotuMPI, zdotcMPI, zaxpy, zcopy, zscal
from .komega_vals_c import get_global_vals_c, initialize_vals_c, cleanup_vals_c
from .komega_vecs_c import get_global_vecs_c, initialize_vecs_c, cleanup_vecs_c


class KomegaCOCG:
    """
    COCG (Conjugate Orthogonal Conjugate Gradient) solver for complex linear systems.

    This class implements the COCG algorithm with shifted equations
    for solving multiple frequency points simultaneously.
    """

    def __init__(self):
        """Initialize the COCG solver."""
        self.params = get_global_params()
        self.math = get_global_math()
        self.vals_c = get_global_vals_c()
        self.vecs_c = get_global_vecs_c()

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
        Initialize the COCG solver.

        Parameters
        ----------
        ndim : int
            Dimension of the Hamiltonian matrix
        nl : int
            Dimension of the projection space
        nz : int
            Number of frequency points
        z : np.ndarray
            Frequency array (complex)
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
        initialize_vals_c(z, itermax)
        initialize_vecs_c(ndim, nl, nz, itermax)

        # Initialize solution array
        x = np.zeros((nl, nz), dtype=complex)

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
                1.0 + self.vals_c.alpha * (self.vals_c.z[iz] - self.vals_c.z_seed)
            ) * self.vals_c.pi[
                iz
            ] - self.vals_c.alpha * self.vals_c.beta / self.vals_c.alpha_old * (
                self.vals_c.pi_old[iz] - self.vals_c.pi[iz]
            )

            # Update search direction
            self.vecs_c.p[:, iz] = (
                r_l / self.vals_c.pi[iz]
                + (self.vals_c.pi_old[iz] / self.vals_c.pi[iz]) ** 2
                * self.vals_c.beta
                * self.vecs_c.p[:, iz]
            )

            # Update solution
            zaxpy(
                self.vals_c.pi[iz] / pi_new * self.vals_c.alpha,
                self.vecs_c.p[:, iz],
                x[:, iz],
            )

            # Update π values
            self.vals_c.pi_old[iz] = self.vals_c.pi[iz]
            self.vals_c.pi[iz] = pi_new

            # Save for restart if needed
            if self.params.itermax > 0:
                self.vals_c.pi_save[iz, self.params.iter] = pi_new

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
        status[2] = self.vals_c.find_minimum_pi(self.params.lz_conv)

        if abs(self.vals_c.pi[status[2]]) < self.params.almost0:
            status[1] = 3

        if status[2] != self.params.iz_seed:
            # Update seed frequency
            self.params.iz_seed = status[2]
            self.vals_c.z_seed = self.vals_c.z[self.params.iz_seed]

            # Update α and ρ
            self.vals_c.alpha *= (
                self.vals_c.pi_old[self.params.iz_seed]
                / self.vals_c.pi[self.params.iz_seed]
            )
            self.vals_c.rho /= self.vals_c.pi_old[self.params.iz_seed] ** 2

            # Scale vectors
            scale = 1.0 / self.vals_c.pi[self.params.iz_seed]
            zscal(scale, v2)
            self.vals_c.scale_pi_values(scale)

            # Scale old vectors
            scale_old = 1.0 / self.vals_c.pi_old[self.params.iz_seed]
            zscal(scale_old, self.vecs_c.v3)
            self.vals_c.pi_old *= scale_old

    def update(
        self,
        v12: np.ndarray,
        v2: np.ndarray,
        x: np.ndarray,
        r_l: np.ndarray,
        status: List[int],
    ) -> None:
        """
        Update the COCG iteration.

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
        rho_old = self.vals_c.rho
        self.vals_c.rho = zdotuMPI(v2, v2)

        if self.params.iter == 1:
            self.vals_c.beta = 0.0 + 0.0j
        else:
            self.vals_c.beta = self.vals_c.rho / rho_old

        # Update working vector
        v12[:] = self.vals_c.z_seed * v2 - v12

        # Update α
        self.vals_c.alpha_old = self.vals_c.alpha
        alpha_denom = (
            zdotuMPI(v2, v12) - self.vals_c.beta * self.vals_c.rho / self.vals_c.alpha
        )

        if abs(alpha_denom) < self.params.almost0:
            status[1] = 2
        elif abs(self.vals_c.rho) < self.params.almost0:
            status[1] = 4

        self.vals_c.alpha = self.vals_c.rho / alpha_denom

        # Save for restart
        if self.params.itermax > 0:
            self.vals_c.alpha_save[self.params.iter - 1] = self.vals_c.alpha
            self.vals_c.beta_save[self.params.iter - 1] = self.vals_c.beta
            self.vecs_c.save_r_l(r_l, self.params.iter)

        # Update shifted equations
        self.shifted_equation(r_l, x)

        # Update residual vector
        v12[:] = (
            (1.0 + self.vals_c.alpha * self.vals_c.beta / self.vals_c.alpha_old) * v2
            - self.vals_c.alpha * v12
            - self.vals_c.alpha
            * self.vals_c.beta
            / self.vals_c.alpha_old
            * self.vecs_c.v3
        )
        zcopy(v2, self.vecs_c.v3)
        zcopy(v12, v2)

        # Perform seed switching
        self.seed_switch(v2, status)

        # Check convergence
        v12[0] = np.sqrt(np.real(zdotcMPI(v2, v2)))
        self.params.resnorm = np.real(v12[0])

        # Update convergence flags
        for iz in range(self.params.nz):
            if abs(v12[0] / self.vals_c.pi[iz]) < self.params.threshold:
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
        elif status[1] == 4:
            status[0] = -self.params.iter
        else:
            status[0] = self.params.iter
            status[1] = 0

    def get_coefficients(self) -> Tuple[np.ndarray, np.ndarray, complex, np.ndarray]:
        """
        Get saved coefficients for restart.

        Returns
        -------
        tuple
            (alpha_save, beta_save, z_seed, r_l_save)
        """
        if not self.initialized:
            return np.array([]), np.array([]), 0.0, np.array([])

        alpha_save, beta_save, _ = self.vals_c.get_saved_values(self.params.iter)
        r_l_save = self.vecs_c.get_saved_r_l(self.params.iter)

        return alpha_save, beta_save, self.vals_c.z_seed, r_l_save

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

        return self.vecs_c.get_v3()

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

        return self.params.resnorm / np.abs(self.vals_c.pi)

    def finalize(self) -> None:
        """Finalize the COCG solver and clean up memory."""
        if self.initialized:
            cleanup_vals_c()
            cleanup_vecs_c()
            self.initialized = False


# Global COCG solver instance
_global_cocg = KomegaCOCG()


def get_global_cocg() -> KomegaCOCG:
    """
    Get the global COCG solver instance.

    Returns
    -------
    KomegaCOCG
        Global COCG solver instance
    """
    return _global_cocg


# Convenience functions for direct access
def komega_COCG_init(
    ndim: int,
    nl: int,
    nz: int,
    z: np.ndarray,
    itermax: int,
    threshold: float,
    comm: Optional[int] = None,
) -> np.ndarray:
    """Initialize the COCG solver."""
    return _global_cocg.init(ndim, nl, nz, z, itermax, threshold, comm)


def komega_COCG_update(
    v12: np.ndarray, v2: np.ndarray, x: np.ndarray, r_l: np.ndarray, status: List[int]
) -> None:
    """Update the COCG iteration."""
    _global_cocg.update(v12, v2, x, r_l, status)


def komega_COCG_getcoef() -> Tuple[np.ndarray, np.ndarray, complex, np.ndarray]:
    """Get saved coefficients for restart."""
    return _global_cocg.get_coefficients()


def komega_COCG_getvec() -> np.ndarray:
    """Get old residual vector."""
    return _global_cocg.get_vectors()


def komega_COCG_getresidual() -> np.ndarray:
    """Get residual norms for all frequencies."""
    return _global_cocg.get_residual()


def komega_COCG_finalize() -> None:
    """Finalize the COCG solver."""
    _global_cocg.finalize()
