"""
Komega Main Module - Unified interface for the Komega library

This module provides a unified interface for all Komega solvers,
making it easy to use the library from Python applications.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

from typing import List, Optional, Tuple, Union

import numpy as np

from komega_bicg import get_global_bicg
from komega_cg_c import get_global_cg_c
from komega_cg_r import get_global_cg_r
from komega_cocg import get_global_cocg


class KomegaSolver:
    """
    Unified interface for Komega solvers.

    This class provides a unified interface for all Komega solvers,
    making it easy to switch between different algorithms.
    """

    def __init__(self, solver_type: str = "bicg"):
        """
        Initialize the Komega solver.

        Parameters
        ----------
        solver_type : str
            Type of solver to use ('bicg', 'cg_r', 'cg_c', 'cocg')
        """
        self.solver_type = solver_type.lower()
        self.solver = None
        self.initialized = False

        # Initialize the appropriate solver
        if self.solver_type == "bicg":
            self.solver = get_global_bicg()
        elif self.solver_type == "cg_r":
            self.solver = get_global_cg_r()
        elif self.solver_type == "cg_c":
            self.solver = get_global_cg_c()
        elif self.solver_type == "cocg":
            self.solver = get_global_cocg()
        else:
            raise ValueError(f"Unknown solver type: {solver_type}")

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
        Initialize the solver.

        Parameters
        ----------
        ndim : int
            Dimension of the Hamiltonian matrix
        nl : int
            Dimension of the projection space
        nz : int
            Number of frequency points
        z : np.ndarray
            Frequency array
        itermax : int
            Maximum number of iterations
        threshold : float
            Convergence threshold
        comm : int, optional
            MPI communicator (for MPI version)

        Returns
        -------
        np.ndarray
            Initialized solution array
        """
        if self.solver is None:
            raise RuntimeError("Solver not initialized")

        x = self.solver.init(ndim, nl, nz, z, itermax, threshold, comm)
        self.initialized = True
        return x

    def update(self, *args, **kwargs) -> None:
        """
        Update the solver iteration.

        Parameters vary depending on the solver type.
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized")

        if self.solver is None:
            raise RuntimeError("Solver not initialized")

        if self.solver_type == "bicg":
            self.solver.update(*args, **kwargs)
        elif self.solver_type in ["cg_r", "cg_c", "cocg"]:
            self.solver.update(*args, **kwargs)
        else:
            raise RuntimeError(f"Unknown solver type: {self.solver_type}")

    def get_coefficients(
        self,
    ) -> Tuple[np.ndarray, np.ndarray, Union[float, complex], np.ndarray]:
        """
        Get saved coefficients for restart.

        Returns
        -------
        tuple
            (alpha_save, beta_save, z_seed, r_l_save)
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized")

        if self.solver is None:
            raise RuntimeError("Solver not initialized")

        return self.solver.get_coefficients()

    def get_vectors(self) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        Get old residual vectors.

        Returns
        -------
        np.ndarray or tuple
            Old residual vector(s)
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized")

        if self.solver is None:
            raise RuntimeError("Solver not initialized")

        return self.solver.get_vectors()

    def get_residual(self) -> np.ndarray:
        """
        Get residual norms for all frequencies.

        Returns
        -------
        np.ndarray
            Residual norms
        """
        if not self.initialized:
            raise RuntimeError("Solver not initialized")

        if self.solver is None:
            raise RuntimeError("Solver not initialized")

        return self.solver.get_residual()

    def finalize(self) -> None:
        """Finalize the solver and clean up memory."""
        if self.initialized and self.solver is not None:
            self.solver.finalize()
            self.initialized = False


def create_solver(solver_type: str = "bicg") -> KomegaSolver:
    """
    Create a new Komega solver instance.

    Parameters
    ----------
    solver_type : str
        Type of solver to create ('bicg', 'cg_r', 'cg_c', 'cocg')

    Returns
    -------
    KomegaSolver
        New solver instance
    """
    return KomegaSolver(solver_type)


def get_available_solvers() -> List[str]:
    """
    Get list of available solver types.

    Returns
    -------
    List[str]
        List of available solver types
    """
    return ["bicg", "cg_r", "cg_c", "cocg"]


def get_solver_info(solver_type: str) -> dict:
    """
    Get information about a specific solver.

    Parameters
    ----------
    solver_type : str
        Type of solver

    Returns
    -------
    dict
        Solver information
    """
    solver_info = {
        "bicg": {
            "name": "BiCG",
            "description": "Bi-Conjugate Gradient solver for complex linear systems",
            "supports_restart": True,
            "supports_mpi": True,
            "data_type": "complex",
        },
        "cg_r": {
            "name": "CG Real",
            "description": "Conjugate Gradient solver for real linear systems",
            "supports_restart": True,
            "supports_mpi": True,
            "data_type": "real",
        },
        "cg_c": {
            "name": "CG Complex",
            "description": "Conjugate Gradient solver for complex linear systems",
            "supports_restart": True,
            "supports_mpi": True,
            "data_type": "complex",
        },
        "cocg": {
            "name": "COCG",
            "description": "Conjugate Orthogonal Conjugate Gradient solver for complex linear systems",
            "supports_restart": True,
            "supports_mpi": True,
            "data_type": "complex",
        },
    }

    if solver_type not in solver_info:
        raise ValueError(f"Unknown solver type: {solver_type}")

    return solver_info[solver_type]


# Convenience functions for direct access to global solvers
def get_bicg_solver():
    """Get the global BiCG solver instance."""
    return get_global_bicg()


def get_cg_r_solver():
    """Get the global CG real solver instance."""
    return get_global_cg_r()


def get_cg_c_solver():
    """Get the global CG complex solver instance."""
    return get_global_cg_c()


def get_cocg_solver():
    """Get the global COCG solver instance."""
    return get_global_cocg()


# Version information
__version__ = "2.0.0"
__author__ = "Mitsuaki Kawamura"
__email__ = "kawamura@issp.u-tokyo.ac.jp"
__description__ = "Komega - A library for solving linear systems in materials science"
