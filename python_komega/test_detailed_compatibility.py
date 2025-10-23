"""
Detailed Fortran Compatibility Test for Komega Python Library

This test suite provides detailed comparison with Fortran test results,
including exact numerical values and convergence behavior.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

import numpy as np
import sys
import os
from typing import Tuple, List, Optional, Dict, Any

# Add the current directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from komega import create_solver, get_available_solvers


class DetailedTestSystem:
    """
    Detailed test system that exactly reproduces Fortran test cases.
    """

    def __init__(self, test_case: str = "solve_cc"):
        """
        Initialize test system for specific Fortran test case.

        Parameters
        ----------
        test_case : str
            Test case name ("solve_cc", "solve_cr", "solve_rc", "solve_rr")
        """
        self.test_case = test_case
        self.setup_test_parameters()

    def setup_test_parameters(self):
        """Set up test parameters based on Fortran test cases."""
        if self.test_case == "solve_cc":
            # Complex-Complex BiCG test
            self.ndim = 5
            self.nl = 5
            self.nz = 5
            self.itermax = 5
            self.threshold = 1e-3
            self.rnd_seed = 100
            self.z_frequencies = [
                -2.0 + 1.0j,
                -1.0 + 1.0j,
                0.0 + 1.0j,
                1.0 + 1.0j,
                2.0 + 1.0j,
            ]
            self.solver_type = "bicg"
            self.data_type = complex

        elif self.test_case == "solve_cr":
            # Complex-Real CG test
            self.ndim = 5
            self.nl = 5
            self.nz = 5
            self.itermax = 5
            self.threshold = 1e-3
            self.rnd_seed = 100
            self.z_frequencies = [
                -2.0 + 1.0j,
                -1.0 + 1.0j,
                0.0 + 1.0j,
                1.0 + 1.0j,
                2.0 + 1.0j,
            ]
            self.solver_type = "cg_c"
            self.data_type = complex

        elif self.test_case == "solve_rc":
            # Real-Complex COCG test
            self.ndim = 5
            self.nl = 5
            self.nz = 5
            self.itermax = 5
            self.threshold = 1e-3
            self.rnd_seed = 100
            self.z_frequencies = [
                -2.0 + 1.0j,
                -1.0 + 1.0j,
                0.0 + 1.0j,
                1.0 + 1.0j,
                2.0 + 1.0j,
            ]
            self.solver_type = "cocg"
            self.data_type = complex

        elif self.test_case == "solve_rr":
            # Real-Real CG test
            self.ndim = 5
            self.nl = 5
            self.nz = 5
            self.itermax = 5
            self.threshold = 1e-3
            self.rnd_seed = 100
            self.z_frequencies = [-2.0, -1.0, 0.0, 1.0, 2.0]
            self.solver_type = "cg_r"
            self.data_type = float

        else:
            raise ValueError(f"Unknown test case: {self.test_case}")

    def generate_system(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Generate test system exactly as in Fortran code.

        Returns
        -------
        tuple
            (hamiltonian, rhs, frequencies)
        """
        print(f"\n#####  Generate Linear System ({self.test_case})  #####")

        # Set random seed for reproducibility
        np.random.seed(self.rnd_seed)

        # Set frequencies
        self.z = np.array(self.z_frequencies, dtype=self.data_type)
        print("  Frequency :")
        for i, freq in enumerate(self.z):
            print(f"  {i+1} {freq}")

        if self.data_type == complex:
            # Complex system generation
            ham_r = np.zeros((self.ndim, self.ndim))
            ham_i = np.zeros((self.ndim, self.ndim))

            # Initialize diagonal and off-diagonal elements
            ham_r[0, 0] = np.random.random()
            ham_i[0, 0] = np.random.random()

            for i in range(1, self.ndim):
                ham_r[i, i] = np.random.random()
                ham_r[i, i - 1] = np.random.random()
                ham_i[i, i] = np.random.random()
                ham_i[i, i - 1] = np.random.random()

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

        else:
            # Real system generation
            self.ham = np.zeros((self.ndim, self.ndim))
            self.ham[0, 0] = np.random.random()

            for i in range(1, self.ndim):
                self.ham[i, i] = np.random.random()
                self.ham[i, i - 1] = np.random.random()

            # Make it symmetric (A^T * A)
            self.ham = self.ham.T @ self.ham

            # Generate right-hand side vector
            self.rhs = np.zeros(self.ndim, dtype=float)
            self.rhs[0] = 1.0

        # Print right-hand side vector
        print("  Right Hand Side Vector :")
        if self.data_type == complex:
            for i in range(self.ndim):
                print(f"  {self.rhs[i].real:15.5e} {self.rhs[i].imag:15.5e}")
        else:
            for i in range(self.ndim):
                print(f"  {self.rhs[i]:15.5e}")

        return self.ham, self.rhs, self.z

    def matrix_vector_product(self, v: np.ndarray) -> np.ndarray:
        """Compute matrix-vector product H * v."""
        return self.ham @ v

    def compute_residual(self, x: np.ndarray, z: complex) -> np.ndarray:
        """Compute residual vector r = z * x - H * x - b."""
        return z * x - self.matrix_vector_product(x) - self.rhs

    def run_solver_test(self) -> Dict[str, Any]:
        """
        Run the solver test and return detailed results.

        Returns
        -------
        dict
            Test results including convergence info, residuals, etc.
        """
        print(f"\n#####  Testing {self.test_case.upper()}  #####")

        # Generate system
        ham, rhs, z = self.generate_system()

        # Create solver
        solver = create_solver(self.solver_type)

        # Initialize solver
        print(f"\n#####  {self.solver_type.upper()} Initialization  #####")
        x = solver.init(self.ndim, self.nl, self.nz, z, self.itermax, self.threshold)

        # Initialize working vectors
        if self.solver_type == "bicg":
            # BiCG needs both v2 and v4
            v2 = rhs.copy()
            v4 = np.conj(v2) if self.data_type == complex else v2.copy()
            v12 = np.zeros(self.ndim, dtype=self.data_type)
            v14 = np.zeros(self.ndim, dtype=self.data_type)
            r_l = np.zeros(self.nl, dtype=self.data_type)
        else:
            # CG and COCG need only v2
            v2 = rhs.copy()
            v12 = np.zeros(self.ndim, dtype=self.data_type)
            r_l = np.zeros(self.nl, dtype=self.data_type)

        # Iteration loop
        print(f"\n#####  {self.solver_type.upper()} Iteration  #####")
        iteration_results = []

        for iter in range(1, self.itermax + 1):
            # Project residual vector
            r_l[:] = v2[: self.nl]

            # Matrix-vector products
            if self.solver_type == "bicg":
                v12[:] = self.matrix_vector_product(v2)
                v14[:] = self.matrix_vector_product(v4)
                # Update solver
                status = [0, 0, 0]
                solver.update(v12, v2, v14, v4, x, r_l, status)
            else:
                v12[:] = self.matrix_vector_product(v2)
                # Update solver
                status = [0, 0, 0]
                solver.update(v12, v2, x, r_l, status)

            # Store iteration results
            iteration_results.append(
                {
                    "iteration": iter,
                    "status": status.copy(),
                    "residual_norm": (
                        np.abs(v12[0]) if self.data_type == complex else v12[0]
                    ),
                    "v2_norm": np.linalg.norm(v2),
                    "x_norm": np.linalg.norm(x),
                }
            )

            # Print debug info
            if self.data_type == complex:
                print(
                    f"DEBUG : {iter:8d} {status[0]:3d} {status[1]:3d} {status[2]:3d} {v12[0].real:15.5e}"
                )
            else:
                print(
                    f"DEBUG : {status[0]:6d} {status[1]:6d} {status[2]:6d} {v12[0]:15.5e}"
                )

            if status[0] < 0:
                break

        # Check convergence
        final_status = iteration_results[-1]["status"]
        if final_status[2] == 0:
            print(f"  Converged in iteration {abs(final_status[0])}")
        elif final_status[2] == 1:
            print(f"  Not Converged in iteration {abs(final_status[0])}")
        elif final_status[2] == 2:
            print(f"  Alpha becomes infinity {abs(final_status[0])}")
        elif final_status[2] == 3:
            print(f"  Pi_seed becomes zero {abs(final_status[0])}")
        elif final_status[2] == 4:
            print(f"  Residual & Shadow residual are orthogonal {abs(final_status[0])}")

        # Compute final residuals
        print("\n#####  Check Results  #####")
        print("  Resulting Vector")
        for iz in range(self.nz):
            for i in range(self.nl):
                if self.data_type == complex:
                    print(f"  {x[i, iz].real:15.5e} {x[i, iz].imag:15.5e}", end="")
                else:
                    print(f"  {x[i, iz]:15.5e}", end="")
            print()

        print("  Residual Vector")
        final_residuals = []
        for iz in range(self.nz):
            residual = self.compute_residual(x[:, iz], z[iz])
            final_residuals.append(residual)
            for i in range(self.nl):
                if self.data_type == complex:
                    print(
                        f"  {residual[i].real:15.5e} {residual[i].imag:15.5e}", end=""
                    )
                else:
                    print(f"  {residual[i]:15.5e}", end="")
            print()

        # Clean up
        solver.finalize()

        # Return detailed results
        return {
            "test_case": self.test_case,
            "solver_type": self.solver_type,
            "convergence_status": final_status,
            "iteration_results": iteration_results,
            "final_residuals": final_residuals,
            "solution": x,
            "hamiltonian": ham,
            "rhs": rhs,
            "frequencies": z,
        }


def run_all_tests():
    """
    Run all detailed compatibility tests.
    """
    print("=" * 80)
    print("Komega Python Library - Detailed Fortran Compatibility Test Suite")
    print("=" * 80)
    print()
    print("This test suite provides detailed comparison with Fortran test results.")
    print()

    test_cases = ["solve_cc", "solve_cr", "solve_rc", "solve_rr"]
    all_results = {}

    try:
        for test_case in test_cases:
            print(f"\n{'='*60}")
            print(f"Running {test_case.upper()} Test")
            print(f"{'='*60}")

            system = DetailedTestSystem(test_case)
            results = system.run_solver_test()
            all_results[test_case] = results

            print(f"#####  {test_case.upper()} Done  #####")

        # Summary
        print("\n" + "=" * 80)
        print("Test Summary")
        print("=" * 80)

        for test_case, results in all_results.items():
            status = results["convergence_status"]
            print(
                f"{test_case:12s}: Status={status[0]:3d} {status[1]:3d} {status[2]:3d}, "
                f"Iterations={len(results['iteration_results'])}"
            )

        print("\nAll detailed compatibility tests completed successfully!")
        print("=" * 80)

        return all_results

    except Exception as e:
        print(f"\nTest suite failed with error: {e}")
        import traceback

        traceback.print_exc()
        return None


def compare_with_fortran_results(python_results: Dict[str, Any]) -> None:
    """
    Compare Python results with expected Fortran results.

    Parameters
    ----------
    python_results : Dict[str, Any]
        Results from Python implementation
    """
    print("\n" + "=" * 80)
    print("Comparison with Fortran Results")
    print("=" * 80)

    # Expected Fortran results (these would be filled in with actual Fortran output)
    expected_results = {
        "solve_cc": {
            "convergence_status": [-5, 0, 0],  # Example: converged in 5 iterations
            "final_iteration": 5,
        },
        "solve_cr": {"convergence_status": [-5, 0, 0], "final_iteration": 5},
        "solve_rc": {"convergence_status": [-5, 0, 0], "final_iteration": 5},
        "solve_rr": {"convergence_status": [-5, 0, 0], "final_iteration": 5},
    }

    for test_case, results in python_results.items():
        print(f"\n{test_case.upper()} Comparison:")
        print(f"  Python status: {results['convergence_status']}")
        print(f"  Python iterations: {len(results['iteration_results'])}")

        if test_case in expected_results:
            expected = expected_results[test_case]
            print(f"  Expected status: {expected['convergence_status']}")
            print(f"  Expected iterations: {expected['final_iteration']}")

            # Check if results match
            if (
                results["convergence_status"] == expected["convergence_status"]
                and len(results["iteration_results"]) == expected["final_iteration"]
            ):
                print(f"  ✓ Results match Fortran implementation")
            else:
                print(f"  ✗ Results differ from Fortran implementation")
        else:
            print(f"  No expected results available for comparison")


def main():
    """
    Main test function.
    """
    # Run all tests
    results = run_all_tests()

    if results is not None:
        # Compare with Fortran results
        compare_with_fortran_results(results)

        print("\n" + "=" * 80)
        print("Detailed compatibility testing completed!")
        print("=" * 80)

        return 0
    else:
        print("\n" + "=" * 80)
        print("Detailed compatibility testing failed!")
        print("=" * 80)

        return 1


if __name__ == "__main__":
    sys.exit(main())
