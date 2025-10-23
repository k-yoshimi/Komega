"""
Komega Values Complex Module - Complex-valued storage for BiCG/COCG algorithms

This module corresponds to the Fortran komega_vals_c module and provides
storage for complex-valued variables used in BiCG and COCG algorithms.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

from typing import List, Optional

import numpy as np

from .komega_parameter import get_global_params


class KomegaValsC:
    """
    Complex-valued storage for BiCG and COCG algorithms.

    This class manages complex-valued variables used in the BiCG and COCG
    algorithms for solving linear systems with complex coefficients.
    """

    def __init__(self):
        """Initialize complex-valued storage."""
        self.params = get_global_params()

        # Scalar values
        self.z_seed: complex = 0.0 + 0.0j
        self.rho: complex = 0.0 + 0.0j
        self.alpha: complex = 0.0 + 0.0j
        self.alpha_old: complex = 0.0 + 0.0j
        self.beta: complex = 0.0 + 0.0j

        # Arrays
        self.z: Optional[np.ndarray] = None  # Frequencies
        self.pi: Optional[np.ndarray] = None  # π values
        self.pi_old: Optional[np.ndarray] = None  # π values at previous step
        self.pi_save: Optional[np.ndarray] = None  # π values saved for restart
        self.alpha_save: Optional[np.ndarray] = None  # α values saved for restart
        self.beta_save: Optional[np.ndarray] = None  # β values saved for restart

    def initialize(self, z: np.ndarray, itermax: int) -> None:
        """
        Initialize storage arrays.

        Parameters
        ----------
        z : np.ndarray
            Frequency array (complex)
        itermax : int
            Maximum number of iterations
        """
        self.z = z.copy()
        self.nz = len(z)

        # Initialize π arrays
        self.pi = np.ones(self.nz, dtype=complex)
        self.pi_old = np.ones(self.nz, dtype=complex)

        # Initialize scalar values
        self.z_seed = z[0]
        self.rho = 1.0 + 0.0j
        self.alpha = 1.0 + 0.0j
        self.alpha_old = 1.0 + 0.0j
        self.beta = 0.0 + 0.0j

        # Initialize restart arrays if needed
        if itermax > 0:
            self.pi_save = np.ones((self.nz, itermax + 1), dtype=complex)
            self.alpha_save = np.zeros(itermax, dtype=complex)
            self.beta_save = np.zeros(itermax, dtype=complex)

    def update_alpha_beta(self, alpha: complex, beta: complex) -> None:
        """
        Update α and β values.

        Parameters
        ----------
        alpha : complex
            New α value
        beta : complex
            New β value
        """
        self.alpha_old = self.alpha
        self.alpha = alpha
        self.beta = beta

    def save_iteration(
        self, iter_count: int, alpha: complex, beta: complex, pi_values: np.ndarray
    ) -> None:
        """
        Save iteration values for restart.

        Parameters
        ----------
        iter_count : int
            Current iteration number
        alpha : complex
            α value for this iteration
        beta : complex
            β value for this iteration
        pi_values : np.ndarray
            π values for this iteration
        """
        if self.alpha_save is not None and iter_count <= len(self.alpha_save):
            self.alpha_save[iter_count - 1] = alpha
        if self.beta_save is not None and iter_count <= len(self.beta_save):
            self.beta_save[iter_count - 1] = beta
        if self.pi_save is not None and iter_count <= self.pi_save.shape[1]:
            self.pi_save[:, iter_count] = pi_values

    def get_saved_values(self, iter_count: int) -> tuple:
        """
        Get saved values for restart.

        Parameters
        ----------
        iter_count : int
            Number of iterations to retrieve

        Returns
        -------
        tuple
            (alpha_save, beta_save, pi_save) arrays
        """
        if iter_count == 0:
            return np.array([]), np.array([]), np.array([])

        alpha_save = (
            self.alpha_save[:iter_count]
            if self.alpha_save is not None
            else np.array([])
        )
        beta_save = (
            self.beta_save[:iter_count] if self.beta_save is not None else np.array([])
        )
        pi_save = (
            self.pi_save[:, :iter_count] if self.pi_save is not None else np.array([])
        )

        return alpha_save, beta_save, pi_save

    def update_pi(self, z: np.ndarray) -> None:
        """
        Update π values for shifted equations.

        Parameters
        ----------
        z : np.ndarray
            Current frequency values
        """
        if self.pi is None or self.pi_old is None:
            return

        # Update π values using the shifted equation formula
        pi_new = (
            1.0 + self.alpha * (z - self.z_seed)
        ) * self.pi - self.alpha * self.beta / self.alpha_old * (self.pi_old - self.pi)

        self.pi_old = self.pi.copy()
        self.pi = pi_new

    def set_seed_frequency(self, iz: int) -> None:
        """
        Set the seed frequency index.

        Parameters
        ----------
        iz : int
            Index of the seed frequency
        """
        if self.z is not None and 0 <= iz < len(self.z):
            self.z_seed = self.z[iz]

    def scale_pi_values(self, scale_factor: complex) -> None:
        """
        Scale π values by a factor.

        Parameters
        ----------
        scale_factor : complex
            Scaling factor
        """
        if self.pi is not None:
            self.pi *= scale_factor
        if self.pi_old is not None:
            self.pi_old *= scale_factor

    def get_pi_values(self) -> np.ndarray:
        """
        Get current π values.

        Returns
        -------
        np.ndarray
            Current π values
        """
        return self.pi.copy() if self.pi is not None else np.array([])

    def get_pi_old_values(self) -> np.ndarray:
        """
        Get previous π values.

        Returns
        -------
        np.ndarray
            Previous π values
        """
        return self.pi_old.copy() if self.pi_old is not None else np.array([])

    def find_minimum_pi(self, lz_conv: List[bool]) -> int:
        """
        Find the index of the minimum |π| value among unconverged frequencies.

        Parameters
        ----------
        lz_conv : List[bool]
            Convergence flags for each frequency

        Returns
        -------
        int
            Index of minimum |π| value
        """
        if self.pi is None:
            return 0

        # Find unconverged frequencies
        unconverged_indices = [i for i, conv in enumerate(lz_conv) if not conv]

        if not unconverged_indices:
            return 0

        # Find minimum |π| among unconverged frequencies
        pi_abs = np.abs(self.pi[unconverged_indices])
        min_idx = np.argmin(pi_abs)

        return unconverged_indices[min_idx]

    def cleanup(self) -> None:
        """Clean up allocated arrays."""
        self.z = None
        self.pi = None
        self.pi_old = None
        self.pi_save = None
        self.alpha_save = None
        self.beta_save = None


# Global complex values instance
_global_vals_c = KomegaValsC()


def get_global_vals_c() -> KomegaValsC:
    """
    Get the global complex values instance.

    Returns
    -------
    KomegaValsC
        Global complex values instance
    """
    return _global_vals_c


def initialize_vals_c(z: np.ndarray, itermax: int) -> None:
    """
    Initialize global complex values.

    Parameters
    ----------
    z : np.ndarray
        Frequency array (complex)
    itermax : int
        Maximum number of iterations
    """
    _global_vals_c.initialize(z, itermax)


def cleanup_vals_c() -> None:
    """Clean up global complex values."""
    _global_vals_c.cleanup()
