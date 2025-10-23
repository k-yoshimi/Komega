"""
Komega Vectors Complex Module - Complex vector storage for BiCG/COCG algorithms

This module corresponds to the Fortran komega_vecs_c module and provides
storage for complex vectors used in BiCG and COCG algorithms.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

from typing import Optional

import numpy as np

from komega_parameter import get_global_params


class KomegaVecsC:
    """
    Complex vector storage for BiCG and COCG algorithms.

    This class manages complex vectors used in the BiCG and COCG
    algorithms for solving linear systems with complex coefficients.
    """

    def __init__(self):
        """Initialize complex vector storage."""
        self.params = get_global_params()

        # Working vectors
        self.v3: Optional[np.ndarray] = None  # Working vector 3
        self.v5: Optional[np.ndarray] = None  # Working vector 5
        self.p: Optional[np.ndarray] = None  # Search direction vectors
        self.r_l_save: Optional[np.ndarray] = None  # Saved residual vectors

    def initialize(self, ndim: int, nl: int, nz: int, itermax: int) -> None:
        """
        Initialize vector storage arrays.

        Parameters
        ----------
        ndim : int
            Dimension of the Hamiltonian matrix
        nl : int
            Dimension of the projection space
        nz : int
            Number of frequency points
        itermax : int
            Maximum number of iterations
        """
        # Initialize working vectors
        self.v3 = np.zeros(ndim, dtype=complex)
        self.v5 = np.zeros(ndim, dtype=complex)
        self.p = np.zeros((nl, nz), dtype=complex)

        # Initialize restart vectors if needed
        if itermax > 0:
            self.r_l_save = np.zeros((nl, itermax), dtype=complex)

    def set_v3(self, v: np.ndarray) -> None:
        """
        Set the v3 working vector.

        Parameters
        ----------
        v : np.ndarray
            Vector to set
        """
        if self.v3 is not None and len(v) == len(self.v3):
            self.v3[:] = v[:]

    def get_v3(self) -> np.ndarray:
        """
        Get the v3 working vector.

        Returns
        -------
        np.ndarray
            Copy of the v3 vector
        """
        return self.v3.copy() if self.v3 is not None else np.array([])

    def set_v5(self, v: np.ndarray) -> None:
        """
        Set the v5 working vector.

        Parameters
        ----------
        v : np.ndarray
            Vector to set
        """
        if self.v5 is not None and len(v) == len(self.v5):
            self.v5[:] = v[:]

    def get_v5(self) -> np.ndarray:
        """
        Get the v5 working vector.

        Returns
        -------
        np.ndarray
            Copy of the v5 vector
        """
        return self.v5.copy() if self.v5 is not None else np.array([])

    def set_p(self, p_values: np.ndarray, iz: int) -> None:
        """
        Set the search direction vector for frequency iz.

        Parameters
        ----------
        p_values : np.ndarray
            Search direction vector
        iz : int
            Frequency index
        """
        if self.p is not None and 0 <= iz < self.p.shape[1]:
            self.p[:, iz] = p_values

    def get_p(self, iz: int) -> np.ndarray:
        """
        Get the search direction vector for frequency iz.

        Parameters
        ----------
        iz : int
            Frequency index

        Returns
        -------
        np.ndarray
            Search direction vector
        """
        if self.p is not None and 0 <= iz < self.p.shape[1]:
            return self.p[:, iz].copy()
        return np.array([])

    def save_r_l(self, r_l: np.ndarray, iter_count: int) -> None:
        """
        Save residual vector for restart.

        Parameters
        ----------
        r_l : np.ndarray
            Residual vector
        iter_count : int
            Iteration number
        """
        if self.r_l_save is not None and iter_count <= self.r_l_save.shape[1]:
            self.r_l_save[:, iter_count - 1] = r_l

    def get_saved_r_l(self, iter_count: int) -> np.ndarray:
        """
        Get saved residual vectors.

        Parameters
        ----------
        iter_count : int
            Number of iterations to retrieve

        Returns
        -------
        np.ndarray
            Saved residual vectors
        """
        if self.r_l_save is not None and iter_count > 0:
            return self.r_l_save[:, :iter_count].copy()
        return np.array([])

    def scale_v3(self, scale_factor: complex) -> None:
        """
        Scale the v3 vector by a factor.

        Parameters
        ----------
        scale_factor : complex
            Scaling factor
        """
        if self.v3 is not None:
            self.v3 *= scale_factor

    def scale_v5(self, scale_factor: complex) -> None:
        """
        Scale the v5 vector by a factor.

        Parameters
        ----------
        scale_factor : complex
            Scaling factor
        """
        if self.v5 is not None:
            self.v5 *= scale_factor

    def scale_p(self, scale_factor: complex, iz: int) -> None:
        """
        Scale the search direction vector for frequency iz.

        Parameters
        ----------
        scale_factor : complex
            Scaling factor
        iz : int
            Frequency index
        """
        if self.p is not None and 0 <= iz < self.p.shape[1]:
            self.p[:, iz] *= scale_factor

    def scale_all_p(self, scale_factor: complex) -> None:
        """
        Scale all search direction vectors.

        Parameters
        ----------
        scale_factor : complex
            Scaling factor
        """
        if self.p is not None:
            self.p *= scale_factor

    def cleanup(self) -> None:
        """Clean up allocated arrays."""
        self.v3 = None
        self.v5 = None
        self.p = None
        self.r_l_save = None


# Global complex vectors instance
_global_vecs_c = KomegaVecsC()


def get_global_vecs_c() -> KomegaVecsC:
    """
    Get the global complex vectors instance.

    Returns
    -------
    KomegaVecsC
        Global complex vectors instance
    """
    return _global_vecs_c


def initialize_vecs_c(ndim: int, nl: int, nz: int, itermax: int) -> None:
    """
    Initialize global complex vectors.

    Parameters
    ----------
    ndim : int
        Dimension of the Hamiltonian matrix
    nl : int
        Dimension of the projection space
    nz : int
        Number of frequency points
    itermax : int
        Maximum number of iterations
    """
    _global_vecs_c.initialize(ndim, nl, nz, itermax)


def cleanup_vecs_c() -> None:
    """Clean up global complex vectors."""
    _global_vecs_c.cleanup()
