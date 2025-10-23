"""
Komega Vectors Real Module - Real vector storage for CG algorithms

This module corresponds to the Fortran komega_vecs_r module and provides
storage for real vectors used in CG algorithms.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

import numpy as np
from typing import Optional
from komega_parameter import get_global_params


class KomegaVecsR:
    """
    Real vector storage for CG algorithms.
    
    This class manages real vectors used in the CG algorithm
    for solving linear systems with real coefficients.
    """
    
    def __init__(self):
        """Initialize real vector storage."""
        self.params = get_global_params()
        
        # Working vectors
        self.v3: Optional[np.ndarray] = None        # Working vector 3
        self.p: Optional[np.ndarray] = None         # Search direction vectors
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
        self.v3 = np.zeros(ndim, dtype=float)
        self.p = np.zeros((nl, nz), dtype=float)
        
        # Initialize restart vectors if needed
        if itermax > 0:
            self.r_l_save = np.zeros((nl, itermax), dtype=float)
    
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
    
    def scale_v3(self, scale_factor: float) -> None:
        """
        Scale the v3 vector by a factor.
        
        Parameters
        ----------
        scale_factor : float
            Scaling factor
        """
        if self.v3 is not None:
            self.v3 *= scale_factor
    
    def scale_p(self, scale_factor: float, iz: int) -> None:
        """
        Scale the search direction vector for frequency iz.
        
        Parameters
        ----------
        scale_factor : float
            Scaling factor
        iz : int
            Frequency index
        """
        if self.p is not None and 0 <= iz < self.p.shape[1]:
            self.p[:, iz] *= scale_factor
    
    def scale_all_p(self, scale_factor: float) -> None:
        """
        Scale all search direction vectors.
        
        Parameters
        ----------
        scale_factor : float
            Scaling factor
        """
        if self.p is not None:
            self.p *= scale_factor
    
    def cleanup(self) -> None:
        """Clean up allocated arrays."""
        self.v3 = None
        self.p = None
        self.r_l_save = None


# Global real vectors instance
_global_vecs_r = KomegaVecsR()


def get_global_vecs_r() -> KomegaVecsR:
    """
    Get the global real vectors instance.
    
    Returns
    -------
    KomegaVecsR
        Global real vectors instance
    """
    return _global_vecs_r


def initialize_vecs_r(ndim: int, nl: int, nz: int, itermax: int) -> None:
    """
    Initialize global real vectors.
    
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
    _global_vecs_r.initialize(ndim, nl, nz, itermax)


def cleanup_vecs_r() -> None:
    """Clean up global real vectors."""
    _global_vecs_r.cleanup()
