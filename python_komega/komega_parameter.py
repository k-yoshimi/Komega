"""
Komega Parameter Module - Global parameters and constants

This module corresponds to the Fortran komega_parameter module and contains
global parameters used throughout the Komega library.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

import numpy as np
from typing import Optional, List


class KomegaParameter:
    """
    Global parameters for the Komega library.
    
    This class manages global parameters that are used across different
    solver modules in the Komega library.
    """
    
    def __init__(self):
        """Initialize global parameters with default values."""
        # Mathematical constants
        self.almost0 = 1e-50
        
        # MPI and communication parameters
        self.comm: Optional[int] = None
        self.lmpi: bool = False
        
        # Problem dimensions
        self.ndim: int = 0  # Dimension of Hamiltonian
        self.nl: int = 0    # Dimension of projection
        self.nz: int = 0    # Number of frequencies (shifts)
        
        # Iteration parameters
        self.itermax: int = 0  # Maximum number of iterations
        self.iter: int = 0     # Current iteration counter
        self.iz_seed: int = 0  # Index of frequency seed
        
        # Convergence parameters
        self.threshold: float = 1e-6  # Convergence threshold
        self.resnorm: float = 0.0     # Residual norm
        
        # Convergence flags for each frequency
        self.lz_conv: Optional[List[bool]] = None
        
    def initialize(self, ndim: int, nl: int, nz: int, itermax: int, 
                   threshold: float, comm: Optional[int] = None) -> None:
        """
        Initialize parameters for a new problem.
        
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
        threshold : float
            Convergence threshold
        comm : int, optional
            MPI communicator (for MPI version)
        """
        self.ndim = ndim
        self.nl = nl
        self.nz = nz
        self.itermax = itermax
        self.threshold = threshold
        self.iter = 0
        self.iz_seed = 1
        self.resnorm = 0.0
        
        # Initialize convergence flags
        self.lz_conv = [False] * nz
        
        # MPI setup
        if comm is not None:
            self.comm = comm
            self.lmpi = True
        else:
            self.comm = None
            self.lmpi = False
    
    def reset_iteration(self) -> None:
        """Reset iteration counter."""
        self.iter = 0
    
    def increment_iteration(self) -> None:
        """Increment iteration counter."""
        self.iter += 1
    
    def set_seed_frequency(self, iz: int) -> None:
        """
        Set the seed frequency index.
        
        Parameters
        ----------
        iz : int
            Index of the seed frequency
        """
        self.iz_seed = iz
    
    def check_convergence(self, residual_norms: np.ndarray) -> bool:
        """
        Check convergence for all frequencies.
        
        Parameters
        ----------
        residual_norms : np.ndarray
            Residual norms for each frequency
            
        Returns
        -------
        bool
            True if all frequencies have converged
        """
        if self.lz_conv is None:
            return False
            
        # Update convergence flags
        for i in range(self.nz):
            if not self.lz_conv[i] and residual_norms[i] < self.threshold:
                self.lz_conv[i] = True
        
        # Check if all frequencies have converged
        return all(self.lz_conv)
    
    def get_converged_frequencies(self) -> List[int]:
        """
        Get list of converged frequency indices.
        
        Returns
        -------
        List[int]
            List of converged frequency indices
        """
        if self.lz_conv is None:
            return []
        return [i for i, converged in enumerate(self.lz_conv) if converged]
    
    def get_unconverged_frequencies(self) -> List[int]:
        """
        Get list of unconverged frequency indices.
        
        Returns
        -------
        List[int]
            List of unconverged frequency indices
        """
        if self.lz_conv is None:
            return list(range(self.nz))
        return [i for i, converged in enumerate(self.lz_conv) if not converged]


# Global parameter instance
_global_params = KomegaParameter()


def get_global_params() -> KomegaParameter:
    """
    Get the global parameter instance.
    
    Returns
    -------
    KomegaParameter
        Global parameter instance
    """
    return _global_params


def initialize_global_params(ndim: int, nl: int, nz: int, itermax: int, 
                           threshold: float, comm: Optional[int] = None) -> None:
    """
    Initialize global parameters.
    
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
    threshold : float
        Convergence threshold
    comm : int, optional
        MPI communicator (for MPI version)
    """
    _global_params.initialize(ndim, nl, nz, itermax, threshold, comm)


def reset_global_params() -> None:
    """Reset global parameters to default values."""
    global _global_params
    _global_params = KomegaParameter()
