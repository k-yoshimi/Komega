"""
Komega Math Module - Mathematical operations and BLAS/LAPACK wrappers

This module corresponds to the Fortran komega_math module and provides
mathematical operations with MPI support for distributed computing.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

from typing import Optional, Union

import numpy as np

from komega_parameter import get_global_params


class KomegaMath:
    """
    Mathematical operations for the Komega library.

    This class provides BLAS/LAPACK-like operations with MPI support
    for distributed computing environments.
    """

    def __init__(self):
        """Initialize the math operations class."""
        self.params = get_global_params()

    def ddot(self, x: np.ndarray, y: np.ndarray) -> float:
        """
        Compute dot product of two real vectors.

        Parameters
        ----------
        x : np.ndarray
            First vector
        y : np.ndarray
            Second vector

        Returns
        -------
        float
            Dot product of x and y
        """
        return np.dot(x, y)

    def zdotc(self, x: np.ndarray, y: np.ndarray) -> complex:
        """
        Compute conjugate dot product of two complex vectors.

        Parameters
        ----------
        x : np.ndarray
            First vector
        y : np.ndarray
            Second vector

        Returns
        -------
        complex
            Conjugate dot product of x and y
        """
        return np.dot(np.conj(x), y)

    def zdotu(self, x: np.ndarray, y: np.ndarray) -> complex:
        """
        Compute dot product of two complex vectors.

        Parameters
        ----------
        x : np.ndarray
            First vector
        y : np.ndarray
            Second vector

        Returns
        -------
        complex
            Dot product of x and y
        """
        return np.dot(x, y)

    def dscal(self, alpha: float, x: np.ndarray) -> None:
        """
        Scale a real vector by a scalar.

        Parameters
        ----------
        alpha : float
            Scaling factor
        x : np.ndarray
            Vector to scale (modified in place)
        """
        x *= alpha

    def zscal(self, alpha: complex, x: np.ndarray) -> None:
        """
        Scale a complex vector by a scalar.

        Parameters
        ----------
        alpha : complex
            Scaling factor
        x : np.ndarray
            Vector to scale (modified in place)
        """
        x *= alpha

    def dcopy(self, x: np.ndarray, y: np.ndarray) -> None:
        """
        Copy a real vector.

        Parameters
        ----------
        x : np.ndarray
            Source vector
        y : np.ndarray
            Destination vector (modified in place)
        """
        y[:] = x[:]

    def zcopy(self, x: np.ndarray, y: np.ndarray) -> None:
        """
        Copy a complex vector.

        Parameters
        ----------
        x : np.ndarray
            Source vector
        y : np.ndarray
            Destination vector (modified in place)
        """
        y[:] = x[:]

    def daxpy(self, alpha: float, x: np.ndarray, y: np.ndarray) -> None:
        """
        Compute y = alpha * x + y for real vectors.

        Parameters
        ----------
        alpha : float
            Scaling factor
        x : np.ndarray
            First vector
        y : np.ndarray
            Second vector (modified in place)
        """
        y += alpha * x

    def zaxpy(self, alpha: complex, x: np.ndarray, y: np.ndarray) -> None:
        """
        Compute y = alpha * x + y for complex vectors.

        Parameters
        ----------
        alpha : complex
            Scaling factor
        x : np.ndarray
            First vector
        y : np.ndarray
            Second vector (modified in place)
        """
        y += alpha * x

    def ddotMPI(self, x: np.ndarray, y: np.ndarray) -> float:
        """
        Compute dot product with MPI reduction.

        Parameters
        ----------
        x : np.ndarray
            First vector
        y : np.ndarray
            Second vector

        Returns
        -------
        float
            Dot product with MPI reduction
        """
        result = self.ddot(x, y)

        # In MPI version, this would perform MPI_Allreduce
        if self.params.lmpi:
            # For now, just return the local result
            # In a real MPI implementation, this would call MPI_Allreduce
            pass

        return result

    def zdotcMPI(self, x: np.ndarray, y: np.ndarray) -> complex:
        """
        Compute conjugate dot product with MPI reduction.

        Parameters
        ----------
        x : np.ndarray
            First vector
        y : np.ndarray
            Second vector

        Returns
        -------
        complex
            Conjugate dot product with MPI reduction
        """
        result = self.zdotc(x, y)

        # In MPI version, this would perform MPI_Allreduce
        if self.params.lmpi:
            # For now, just return the local result
            # In a real MPI implementation, this would call MPI_Allreduce
            pass

        return result

    def zdotuMPI(self, x: np.ndarray, y: np.ndarray) -> complex:
        """
        Compute dot product with MPI reduction.

        Parameters
        ----------
        x : np.ndarray
            First vector
        y : np.ndarray
            Second vector

        Returns
        -------
        complex
            Dot product with MPI reduction
        """
        result = self.zdotu(x, y)

        # In MPI version, this would perform MPI_Allreduce
        if self.params.lmpi:
            # For now, just return the local result
            # In a real MPI implementation, this would call MPI_Allreduce
            pass

        return result

    def dabsmax(self, array: np.ndarray) -> float:
        """
        Find maximum absolute value with MPI reduction.

        Parameters
        ----------
        array : np.ndarray
            Input array

        Returns
        -------
        float
            Maximum absolute value
        """
        result: float = np.max(np.abs(array))

        # In MPI version, this would perform MPI_Allreduce with MPI_MAX
        if self.params.lmpi:
            # For now, just return the local result
            # In a real MPI implementation, this would call MPI_Allreduce
            pass

        return result

    def norm(self, x: np.ndarray) -> float:
        """
        Compute the Euclidean norm of a vector.

        Parameters
        ----------
        x : np.ndarray
            Input vector

        Returns
        -------
        float
            Euclidean norm
        """
        return np.sqrt(np.real(self.zdotcMPI(x, x)))

    def complex_norm(self, x: np.ndarray) -> float:
        """
        Compute the Euclidean norm of a complex vector.

        Parameters
        ----------
        x : np.ndarray
            Input complex vector

        Returns
        -------
        float
            Euclidean norm
        """
        return np.sqrt(np.real(self.zdotcMPI(x, x)))


# Global math operations instance
_global_math = KomegaMath()


def get_global_math() -> KomegaMath:
    """
    Get the global math operations instance.

    Returns
    -------
    KomegaMath
        Global math operations instance
    """
    return _global_math


# Convenience functions for direct access
def ddot(x: np.ndarray, y: np.ndarray) -> float:
    """Compute dot product of two real vectors."""
    return _global_math.ddot(x, y)


def zdotc(x: np.ndarray, y: np.ndarray) -> complex:
    """Compute conjugate dot product of two complex vectors."""
    return _global_math.zdotc(x, y)


def zdotu(x: np.ndarray, y: np.ndarray) -> complex:
    """Compute dot product of two complex vectors."""
    return _global_math.zdotu(x, y)


def dscal(alpha: float, x: np.ndarray) -> None:
    """Scale a real vector by a scalar."""
    _global_math.dscal(alpha, x)


def zscal(alpha: complex, x: np.ndarray) -> None:
    """Scale a complex vector by a scalar."""
    _global_math.zscal(alpha, x)


def dcopy(x: np.ndarray, y: np.ndarray) -> None:
    """Copy a real vector."""
    _global_math.dcopy(x, y)


def zcopy(x: np.ndarray, y: np.ndarray) -> None:
    """Copy a complex vector."""
    _global_math.zcopy(x, y)


def daxpy(alpha: float, x: np.ndarray, y: np.ndarray) -> None:
    """Compute y = alpha * x + y for real vectors."""
    _global_math.daxpy(alpha, x, y)


def zaxpy(alpha: complex, x: np.ndarray, y: np.ndarray) -> None:
    """Compute y = alpha * x + y for complex vectors."""
    _global_math.zaxpy(alpha, x, y)


def ddotMPI(x: np.ndarray, y: np.ndarray) -> float:
    """Compute dot product with MPI reduction."""
    return _global_math.ddotMPI(x, y)


def zdotcMPI(x: np.ndarray, y: np.ndarray) -> complex:
    """Compute conjugate dot product with MPI reduction."""
    return _global_math.zdotcMPI(x, y)


def zdotuMPI(x: np.ndarray, y: np.ndarray) -> complex:
    """Compute dot product with MPI reduction."""
    return _global_math.zdotuMPI(x, y)


def dabsmax(array: np.ndarray) -> float:
    """Find maximum absolute value with MPI reduction."""
    return _global_math.dabsmax(array)
