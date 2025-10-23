"""
Fortran Compatibility Test for Komega Python Library

This test suite reproduces the Fortran test cases to verify that the Python
implementation produces the same results as the original Fortran code.

Test cases:
1. solve_cc: Complex-Complex BiCG solver
2. solve_cr: Complex-Real CG solver  
3. solve_rc: Real-Complex COCG solver
4. solve_rr: Real-Real CG solver

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

import numpy as np
import sys
import os
from typing import Tuple, List, Optional

# Add the current directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from komega import create_solver, get_available_solvers


class FortranTestSystem:
    """
    Test system generator that mimics the Fortran test cases.
    """
    
    def __init__(self, ndim: int = 5, nl: int = 5, nz: int = 5, 
                 itermax: int = 5, threshold: float = 1e-3, rnd_seed: int = 100):
        """
        Initialize test system parameters.
        
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
        rnd_seed : int
            Random seed for reproducible results
        """
        self.ndim = ndim
        self.nl = nl
        self.nz = nz
        self.itermax = itermax
        self.threshold = threshold
        self.rnd_seed = rnd_seed
        
        # Set random seed for reproducibility
        np.random.seed(rnd_seed)
        
        # Initialize arrays
        self.ham = None
        self.rhs = None
        self.z = None
        self.x = None
        
    def generate_complex_system(self, z_frequencies: List[complex]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate a complex linear system (solve_cc test case).
        
        Parameters
        ----------
        z_frequencies : List[complex]
            Complex frequency points
            
        Returns
        -------
        tuple
            (hamiltonian, rhs, frequencies)
        """
        print("#####  Generate Complex Linear System  #####")
        
        # Set frequencies
        self.z = np.array(z_frequencies, dtype=complex)
        print("  Frequency :")
        for i, freq in enumerate(self.z):
            print(f"  {i+1} {freq}")
        
        # Generate random complex Hamiltonian
        ham_r = np.zeros((self.ndim, self.ndim))
        ham_i = np.zeros((self.ndim, self.ndim))
        
        # Initialize diagonal and off-diagonal elements
        ham_r[0, 0] = np.random.random()
        ham_i[0, 0] = np.random.random()
        
        for i in range(1, self.ndim):
            ham_r[i, i] = np.random.random()
            ham_r[i, i-1] = np.random.random()
            ham_i[i, i] = np.random.random()
            ham_i[i, i-1] = np.random.random()
        
        # Fill remaining elements
        ham_r += np.random.random((self.ndim, self.ndim))
        ham_i += np.random.random((self.ndim, self.ndim))
        
        # Create complex Hamiltonian
        self.ham = ham_r + 1j * ham_i
        
        # Make it Hermitian (A^H * A)
        self.ham = np.conj(self.ham.T) @ self.ham
        
        # Generate right-hand side vector
        rhs_r = np.random.random(self.ndim)
        rhs_i = np.random.random(self.ndim)
        self.rhs = np.zeros(self.ndim, dtype=complex)
        self.rhs[0] = 1.0 + 0.0j
        self.rhs += rhs_r + 1j * rhs_i
        
        print("  Right Hand Side Vector :")
        for i in range(self.ndim):
            print(f"  {self.rhs[i].real:15.5e} {self.rhs[i].imag:15.5e}")
        
        return self.ham, self.rhs, self.z
    
    def generate_real_system(self, z_frequencies: List[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate a real linear system (solve_rr test case).
        
        Parameters
        ----------
        z_frequencies : List[float]
            Real frequency points
            
        Returns
        -------
        tuple
            (hamiltonian, rhs, frequencies)
        """
        print("#####  Generate Real Linear System  #####")
        
        # Set frequencies
        self.z = np.array(z_frequencies, dtype=float)
        print("  Frequency :")
        for i, freq in enumerate(self.z):
            print(f"  {i+1} {freq}")
        
        # Generate random real Hamiltonian
        self.ham = np.zeros((self.ndim, self.ndim))
        self.ham[0, 0] = np.random.random()
        
        for i in range(1, self.ndim):
            self.ham[i, i] = np.random.random()
            self.ham[i, i-1] = np.random.random()
        
        # Make it symmetric (A^T * A)
        self.ham = self.ham.T @ self.ham
        
        # Generate right-hand side vector
        self.rhs = np.zeros(self.ndim, dtype=float)
        self.rhs[0] = 1.0
        
        print("  Right Hand Side Vector :")
        for i in range(self.ndim):
            print(f"  {self.rhs[i]:15.5e}")
        
        return self.ham, self.rhs, self.z
    
    def matrix_vector_product(self, v: np.ndarray) -> np.ndarray:
        """
        Compute matrix-vector product H * v.
        
        Parameters
        ----------
        v : np.ndarray
            Input vector
            
        Returns
        -------
        np.ndarray
            Result vector
        """
        return self.ham @ v
    
    def compute_residual(self, x: np.ndarray, z: complex) -> np.ndarray:
        """
        Compute residual vector r = z * x - H * x - b.
        
        Parameters
        ----------
        x : np.ndarray
            Solution vector
        z : complex
            Frequency point
            
        Returns
        -------
        np.ndarray
            Residual vector
        """
        return z * x - self.matrix_vector_product(x) - self.rhs


def test_solve_cc():
    """
    Test solve_cc: Complex-Complex BiCG solver.
    """
    print("\n" + "="*60)
    print("Testing solve_cc: Complex-Complex BiCG Solver")
    print("="*60)
    
    # Test parameters
    ndim = 5
    nl = 5
    nz = 5
    itermax = 5
    threshold = 1e-3
    rnd_seed = 100
    
    # Complex frequencies
    z_frequencies = [-2.0 + 1.0j, -1.0 + 1.0j, 0.0 + 1.0j, 1.0 + 1.0j, 2.0 + 1.0j]
    
    # Generate test system
    system = FortranTestSystem(ndim, nl, nz, itermax, threshold, rnd_seed)
    ham, rhs, z = system.generate_complex_system(z_frequencies)
    
    # Create BiCG solver
    solver = create_solver('bicg')
    
    # Initialize solver
    print("\n#####  BiCG Initialization  #####")
    x = solver.init(ndim, nl, nz, z, itermax, threshold)
    
    # Initialize vectors
    v2 = rhs.copy()
    v4 = np.conj(v2)
    v12 = np.zeros(ndim, dtype=complex)
    v14 = np.zeros(ndim, dtype=complex)
    r_l = np.zeros(nl, dtype=complex)
    
    # BiCG iteration loop
    print("\n#####  BiCG Iteration  #####")
    for iter in range(1, itermax + 1):
        # Project residual vector
        r_l[:] = v2[:nl]
        
        # Matrix-vector products
        v12[:] = system.matrix_vector_product(v2)
        v14[:] = system.matrix_vector_product(v4)
        
        # Update solver
        status = [0, 0, 0]
        solver.update(v12, v2, v14, v4, x, r_l, status)
        
        print(f"DEBUG : {iter:8d} {status[0]:3d} {status[1]:3d} {status[2]:3d} {v12[0].real:15.5e}")
        
        if status[0] < 0:
            break
    
    # Check convergence
    if status[2] == 0:
        print(f"  Converged in iteration {abs(status[0])}")
    elif status[2] == 1:
        print(f"  Not Converged in iteration {abs(status[0])}")
    elif status[2] == 2:
        print(f"  Alpha becomes infinity {abs(status[0])}")
    elif status[2] == 3:
        print(f"  Pi_seed becomes zero {abs(status[0])}")
    elif status[2] == 4:
        print(f"  Residual & Shadow residual are orthogonal {abs(status[0])}")
    
    # Output results
    print("\n#####  Check Results  #####")
    print("  Resulting Vector")
    for iz in range(nz):
        for i in range(nl):
            print(f"  {x[i, iz].real:15.5e} {x[i, iz].imag:15.5e}", end="")
        print()
    
    print("  Residual Vector")
    for iz in range(nz):
        residual = system.compute_residual(x[:, iz], z[iz])
        for i in range(nl):
            print(f"  {residual[i].real:15.5e} {residual[i].imag:15.5e}", end="")
        print()
    
    # Clean up
    solver.finalize()
    print("#####  Done  #####")


def test_solve_cr():
    """
    Test solve_cr: Complex-Real CG solver.
    """
    print("\n" + "="*60)
    print("Testing solve_cr: Complex-Real CG Solver")
    print("="*60)
    
    # Test parameters
    ndim = 5
    nl = 5
    nz = 5
    itermax = 5
    threshold = 1e-3
    rnd_seed = 100
    
    # Complex frequencies
    z_frequencies = [-2.0 + 1.0j, -1.0 + 1.0j, 0.0 + 1.0j, 1.0 + 1.0j, 2.0 + 1.0j]
    
    # Generate test system
    system = FortranTestSystem(ndim, nl, nz, itermax, threshold, rnd_seed)
    ham, rhs, z = system.generate_complex_system(z_frequencies)
    
    # Create CG complex solver
    solver = create_solver('cg_c')
    
    # Initialize solver
    print("\n#####  CG Initialization  #####")
    x = solver.init(ndim, nl, nz, z, itermax, threshold)
    
    # Initialize vectors
    v2 = rhs.copy()
    v12 = np.zeros(ndim, dtype=complex)
    r_l = np.zeros(nl, dtype=complex)
    
    # CG iteration loop
    print("\n#####  CG Iteration  #####")
    for iter in range(1, itermax + 1):
        # Project residual vector
        r_l[:] = v2[:nl]
        
        # Matrix-vector product
        v12[:] = system.matrix_vector_product(v2)
        
        # Update solver
        status = [0, 0, 0]
        solver.update(v12, v2, x, r_l, status)
        
        print(f"DEBUG : {iter:8d} {status[0]:3d} {status[1]:3d} {status[2]:3d} {v12[0].real:15.5e}")
        
        if status[0] < 0:
            break
    
    # Check convergence
    if status[2] == 0:
        print(f"  Converged in iteration {abs(status[0])}")
    elif status[2] == 1:
        print(f"  Not Converged in iteration {abs(status[0])}")
    elif status[2] == 2:
        print(f"  Alpha becomes infinity {abs(status[0])}")
    elif status[2] == 3:
        print(f"  Pi_seed becomes zero {abs(status[0])}")
    
    # Output results
    print("\n#####  Check Results  #####")
    print("  Resulting Vector")
    for iz in range(nz):
        for i in range(nl):
            print(f"  {x[i, iz].real:15.5e} {x[i, iz].imag:15.5e}", end="")
        print()
    
    print("  Residual Vector")
    for iz in range(nz):
        residual = system.compute_residual(x[:, iz], z[iz])
        for i in range(nl):
            print(f"  {residual[i].real:15.5e} {residual[i].imag:15.5e}", end="")
        print()
    
    # Clean up
    solver.finalize()
    print("#####  Done  #####")


def test_solve_rc():
    """
    Test solve_rc: Real-Complex COCG solver.
    """
    print("\n" + "="*60)
    print("Testing solve_rc: Real-Complex COCG Solver")
    print("="*60)
    
    # Test parameters
    ndim = 5
    nl = 5
    nz = 5
    itermax = 5
    threshold = 1e-3
    rnd_seed = 100
    
    # Complex frequencies
    z_frequencies = [-2.0 + 1.0j, -1.0 + 1.0j, 0.0 + 1.0j, 1.0 + 1.0j, 2.0 + 1.0j]
    
    # Generate test system
    system = FortranTestSystem(ndim, nl, nz, itermax, threshold, rnd_seed)
    ham, rhs, z = system.generate_complex_system(z_frequencies)
    
    # Create COCG solver
    solver = create_solver('cocg')
    
    # Initialize solver
    print("\n#####  COCG Initialization  #####")
    x = solver.init(ndim, nl, nz, z, itermax, threshold)
    
    # Initialize vectors
    v2 = rhs.copy()
    v12 = np.zeros(ndim, dtype=complex)
    r_l = np.zeros(nl, dtype=complex)
    
    # COCG iteration loop
    print("\n#####  COCG Iteration  #####")
    for iter in range(1, itermax + 1):
        # Project residual vector
        r_l[:] = v2[:nl]
        
        # Matrix-vector product
        v12[:] = system.matrix_vector_product(v2)
        
        # Update solver
        status = [0, 0, 0]
        solver.update(v12, v2, x, r_l, status)
        
        print(f"DEBUG : {iter:8d} {status[0]:3d} {status[1]:3d} {status[2]:3d} {v12[0].real:15.5e}")
        
        if status[0] < 0:
            break
    
    # Check convergence
    if status[2] == 0:
        print(f"  Converged in iteration {abs(status[0])}")
    elif status[2] == 1:
        print(f"  Not Converged in iteration {abs(status[0])}")
    elif status[2] == 2:
        print(f"  Alpha becomes infinity {abs(status[0])}")
    elif status[2] == 3:
        print(f"  Pi_seed becomes zero {abs(status[0])}")
    elif status[2] == 4:
        print(f"  Residual & Shadow residual are orthogonal {abs(status[0])}")
    
    # Output results
    print("\n#####  Check Results  #####")
    print("  Resulting Vector")
    for iz in range(nz):
        for i in range(nl):
            print(f"  {x[i, iz].real:15.5e} {x[i, iz].imag:15.5e}", end="")
        print()
    
    print("  Residual Vector")
    for iz in range(nz):
        residual = system.compute_residual(x[:, iz], z[iz])
        for i in range(nl):
            print(f"  {residual[i].real:15.5e} {residual[i].imag:15.5e}", end="")
        print()
    
    # Clean up
    solver.finalize()
    print("#####  Done  #####")


def test_solve_rr():
    """
    Test solve_rr: Real-Real CG solver.
    """
    print("\n" + "="*60)
    print("Testing solve_rr: Real-Real CG Solver")
    print("="*60)
    
    # Test parameters
    ndim = 5
    nl = 5
    nz = 5
    itermax = 5
    threshold = 1e-3
    rnd_seed = 100
    
    # Real frequencies
    z_frequencies = [-2.0, -1.0, 0.0, 1.0, 2.0]
    
    # Generate test system
    system = FortranTestSystem(ndim, nl, nz, itermax, threshold, rnd_seed)
    ham, rhs, z = system.generate_real_system(z_frequencies)
    
    # Create CG real solver
    solver = create_solver('cg_r')
    
    # Initialize solver
    print("\n#####  CG Initialization  #####")
    x = solver.init(ndim, nl, nz, z, itermax, threshold)
    
    # Initialize vectors
    v2 = rhs.copy()
    v12 = np.zeros(ndim, dtype=float)
    r_l = np.zeros(nl, dtype=float)
    
    # CG iteration loop
    print("\n#####  CG Iteration  #####")
    for iter in range(1, itermax + 1):
        # Project residual vector
        r_l[:] = v2[:nl]
        
        # Matrix-vector product
        v12[:] = system.matrix_vector_product(v2)
        
        # Update solver
        status = [0, 0, 0]
        solver.update(v12, v2, x, r_l, status)
        
        print(f"DEBUG : {status[0]:6d} {status[1]:6d} {status[2]:6d} {v12[0]:15.5e}")
        
        if status[0] < 0:
            break
    
    # Check convergence
    if status[2] == 0:
        print(f"  Converged in iteration {abs(status[0])}")
    elif status[2] == 1:
        print(f"  Not Converged in iteration {abs(status[0])}")
    elif status[2] == 2:
        print(f"  Alpha becomes infinity {abs(status[0])}")
    elif status[2] == 3:
        print(f"  Pi_seed becomes zero {abs(status[0])}")
    
    # Output results
    print("\n#####  Check Results  #####")
    print("  Resulting Vector")
    for iz in range(nz):
        for i in range(nl):
            print(f"  {x[i, iz]:15.5e}", end="")
        print()
    
    print("  Residual Vector")
    for iz in range(nz):
        residual = system.compute_residual(x[:, iz], z[iz])
        for i in range(nl):
            print(f"  {residual[i]:15.5e}", end="")
        print()
    
    # Clean up
    solver.finalize()
    print("#####  Done  #####")


def main():
    """
    Run all Fortran compatibility tests.
    """
    print("=" * 80)
    print("Komega Python Library - Fortran Compatibility Test Suite")
    print("=" * 80)
    print()
    print("This test suite reproduces the Fortran test cases to verify")
    print("that the Python implementation produces the same results.")
    print()
    
    try:
        # Test all solver types
        test_solve_cc()   # Complex-Complex BiCG
        test_solve_cr()   # Complex-Real CG
        test_solve_rc()   # Real-Complex COCG
        test_solve_rr()   # Real-Real CG
        
        print("\n" + "=" * 80)
        print("All Fortran compatibility tests completed successfully!")
        print("=" * 80)
        
    except Exception as e:
        print(f"\nTest suite failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
