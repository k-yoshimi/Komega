#!/usr/bin/env python3
"""
GitHub Actions Test Suite for Komega Python Library

This script provides comprehensive testing for GitHub Actions CI/CD pipeline.
It includes all necessary tests in a single, maintainable script.

Usage:
    python test_github_actions.py [--test-type TYPE] [--verbose]

Test Types:
    - basic: Basic functionality tests
    - solvers: Solver creation and initialization tests
    - math: Mathematical operations tests
    - storage: Value and vector storage tests
    - compatibility: Fortran compatibility tests
    - all: Run all tests (default)
"""

import argparse
import sys
import traceback
from typing import List, Optional

import numpy as np

# Import Komega modules
try:
    from komega import (
        create_solver,
        get_available_solvers,
        get_global_math,
        get_global_vals_c,
        get_global_vals_r,
        get_global_vecs_c,
        get_global_vecs_r,
        get_solver_info,
    )
    from komega_parameter import get_global_params
except ImportError as e:
    print(f"❌ Import error: {e}")
    sys.exit(1)


class TestResult:
    """Test result container."""

    def __init__(self, name: str, success: bool, message: str = ""):
        self.name = name
        self.success = success
        self.message = message

    def __str__(self):
        status = "✓" if self.success else "✗"
        return f"{status} {self.name}: {self.message}"


class GitHubActionsTestSuite:
    """Comprehensive test suite for GitHub Actions."""

    def __init__(self, verbose: bool = False):
        self.verbose = verbose
        self.results: List[TestResult] = []
        self.failed_tests: List[str] = []

    def log(self, message: str):
        """Log message if verbose mode is enabled."""
        if self.verbose:
            print(f"DEBUG: {message}")

    def add_result(self, name: str, success: bool, message: str = ""):
        """Add a test result."""
        result = TestResult(name, success, message)
        self.results.append(result)
        if not success:
            self.failed_tests.append(name)
        print(result)

    def test_imports(self):
        """Test that all modules can be imported."""
        self.log("Testing imports...")
        try:
            # Test parameter module
            params = get_global_params()
            self.add_result(
                "Import parameters", True, f"Global params: {type(params).__name__}"
            )

            # Test math module
            math_ops = get_global_math()
            self.add_result(
                "Import math", True, f"Math operations: {type(math_ops).__name__}"
            )

            # Test value storage modules
            vals_r = get_global_vals_r()
            vals_c = get_global_vals_c()
            self.add_result(
                "Import value storage",
                True,
                f"Real: {type(vals_r).__name__}, Complex: {type(vals_c).__name__}",
            )

            # Test vector storage modules
            vecs_r = get_global_vecs_r()
            vecs_c = get_global_vecs_c()
            self.add_result(
                "Import vector storage",
                True,
                f"Real: {type(vecs_r).__name__}, Complex: {type(vecs_c).__name__}",
            )

        except Exception as e:
            self.add_result("Import modules", False, f"Error: {e}")

    def test_solver_creation(self):
        """Test solver creation and basic functionality."""
        self.log("Testing solver creation...")
        try:
            solvers = get_available_solvers()
            self.add_result("Get available solvers", True, f"Found: {solvers}")

            for solver_type in solvers:
                try:
                    solver = create_solver(solver_type)
                    self.add_result(
                        f"Create {solver_type} solver",
                        True,
                        f"Type: {type(solver).__name__}",
                    )
                except Exception as e:
                    self.add_result(
                        f"Create {solver_type} solver", False, f"Error: {e}"
                    )

        except Exception as e:
            self.add_result("Solver creation", False, f"Error: {e}")

    def test_solver_initialization(self):
        """Test solver initialization with different parameters."""
        self.log("Testing solver initialization...")

        # Test parameters
        ndim, nl, nz = 5, 5, 3
        itermax, threshold = 100, 1e-6

        # Test real solver
        try:
            solver_r = create_solver("cg_r")
            z_real = np.array([1.0, 2.0, 3.0])
            x_r = solver_r.init(ndim, nl, nz, z_real, itermax, threshold)
            self.add_result("Real CG initialization", True, f"Shape: {x_r.shape}")
            solver_r.finalize()
        except Exception as e:
            self.add_result("Real CG initialization", False, f"Error: {e}")

        # Test complex solver
        try:
            solver_c = create_solver("bicg")
            z_complex = np.array([1.0 + 0.1j, 2.0 + 0.2j, 3.0 + 0.3j])
            x_c = solver_c.init(ndim, nl, nz, z_complex, itermax, threshold)
            self.add_result("Complex BiCG initialization", True, f"Shape: {x_c.shape}")
            solver_c.finalize()
        except Exception as e:
            self.add_result("Complex BiCG initialization", False, f"Error: {e}")

    def test_mathematical_operations(self):
        """Test mathematical operations."""
        self.log("Testing mathematical operations...")
        try:
            math_ops = get_global_math()

            # Test vectors
            n = 5
            x_real = np.random.rand(n)
            y_real = np.random.rand(n)
            x_complex = np.random.rand(n) + 1j * np.random.rand(n)
            y_complex = np.random.rand(n) + 1j * np.random.rand(n)

            # Test dot products
            dot_real = math_ops.ddot(x_real, y_real)
            dotc_complex = math_ops.zdotc(x_complex, y_complex)
            dotu_complex = math_ops.zdotu(x_complex, y_complex)

            self.add_result("Real dot product", True, f"Result: {dot_real:.6f}")
            self.add_result(
                "Complex conjugate dot", True, f"Result: {dotc_complex:.6f}"
            )
            self.add_result("Complex dot product", True, f"Result: {dotu_complex:.6f}")

            # Test scaling
            alpha_real = 2.0
            alpha_complex = 1.0 + 1j

            x_scaled = x_real.copy()
            math_ops.dscal(alpha_real, x_scaled)
            self.add_result("Real scaling", True, "Operation completed")

            x_complex_scaled = x_complex.copy()
            math_ops.zscal(alpha_complex, x_complex_scaled)
            self.add_result("Complex scaling", True, "Operation completed")

        except Exception as e:
            self.add_result("Mathematical operations", False, f"Error: {e}")

    def test_value_storage(self):
        """Test value storage modules."""
        self.log("Testing value storage...")
        try:
            # Test real values
            vals_r = get_global_vals_r()
            z_real = np.array([1.0, 2.0, 3.0])
            vals_r.initialize(z_real, 100)
            pi_shape = vals_r.get_pi_values().shape
            self.add_result("Real values storage", True, f"π shape: {pi_shape}")

            # Test complex values
            vals_c = get_global_vals_c()
            z_complex = np.array([1.0 + 0.1j, 2.0 + 0.2j, 3.0 + 0.3j])
            vals_c.initialize(z_complex, 100)
            pi_shape = vals_c.get_pi_values().shape
            self.add_result("Complex values storage", True, f"π shape: {pi_shape}")

        except Exception as e:
            self.add_result("Value storage", False, f"Error: {e}")

    def test_vector_storage(self):
        """Test vector storage modules."""
        self.log("Testing vector storage...")
        try:
            # Test real vectors
            vecs_r = get_global_vecs_r()
            vecs_r.initialize(10, 5, 3, 100)
            v3_shape = vecs_r.get_v3().shape
            p_shape = vecs_r.get_p(0).shape
            self.add_result(
                "Real vector storage", True, f"v3: {v3_shape}, p: {p_shape}"
            )

            # Test complex vectors
            vecs_c = get_global_vecs_c()
            vecs_c.initialize(10, 5, 3, 100)
            v3_shape = vecs_c.get_v3().shape
            p_shape = vecs_c.get_p(0).shape
            self.add_result(
                "Complex vector storage", True, f"v3: {v3_shape}, p: {p_shape}"
            )

        except Exception as e:
            self.add_result("Vector storage", False, f"Error: {e}")

    def test_solver_information(self):
        """Test solver information retrieval."""
        self.log("Testing solver information...")
        try:
            solvers = get_available_solvers()
            self.add_result("Available solvers", True, f"Count: {len(solvers)}")

            for solver_type in solvers:
                try:
                    info = get_solver_info(solver_type)
                    self.add_result(
                        f"Solver info {solver_type}",
                        True,
                        f"Name: {info['name']}, Type: {info['data_type']}",
                    )
                except Exception as e:
                    self.add_result(f"Solver info {solver_type}", False, f"Error: {e}")

        except Exception as e:
            self.add_result("Solver information", False, f"Error: {e}")

    def test_error_handling(self):
        """Test error handling."""
        self.log("Testing error handling...")
        try:
            # Test invalid solver type
            try:
                solver = create_solver("invalid_solver")
                self.add_result(
                    "Invalid solver type", False, "Should have raised ValueError"
                )
            except ValueError:
                self.add_result(
                    "Invalid solver type", True, "Correctly raised ValueError"
                )

            # Test solver operations without initialization
            try:
                solver = create_solver("bicg")
                solver.update([], [], [], [], [], [], [])
                self.add_result(
                    "Uninitialized solver", False, "Should have raised RuntimeError"
                )
            except RuntimeError:
                self.add_result(
                    "Uninitialized solver", True, "Correctly raised RuntimeError"
                )

        except Exception as e:
            self.add_result("Error handling", False, f"Error: {e}")

    def test_memory_management(self):
        """Test memory management."""
        self.log("Testing memory management...")
        try:
            # Test multiple solver instances
            solvers = []
            for i in range(5):
                solver = create_solver("cg_r")
                z = np.array([1.0, 2.0, 3.0])
                x = solver.init(5, 5, 3, z, 100, 1e-6)
                solvers.append(solver)

            self.add_result(
                "Multiple solver instances", True, f"Created: {len(solvers)}"
            )

            # Clean up
            for solver in solvers:
                solver.finalize()

            self.add_result("Solver cleanup", True, "All solvers finalized")

        except Exception as e:
            self.add_result("Memory management", False, f"Error: {e}")

    def test_numerical_stability(self):
        """Test numerical stability."""
        self.log("Testing numerical stability...")
        try:
            # Test with small numbers
            solver = create_solver("cg_r")
            z = np.array([1e-10, 1e-8, 1e-6])
            x = solver.init(3, 3, 3, z, 100, 1e-12)
            self.add_result("Small numbers", True, "Solver handles small numbers")
            solver.finalize()

            # Test with large numbers
            solver = create_solver("cg_r")
            z = np.array([1e6, 1e8, 1e10])
            x = solver.init(3, 3, 3, z, 100, 1e-6)
            self.add_result("Large numbers", True, "Solver handles large numbers")
            solver.finalize()

        except Exception as e:
            self.add_result("Numerical stability", False, f"Error: {e}")

    def run_tests(self, test_type: str = "all"):
        """Run specified tests."""
        print(f"🚀 Running GitHub Actions Test Suite: {test_type}")
        print("=" * 60)

        if test_type in ["basic", "all"]:
            print("\n📦 Testing Basic Functionality")
            print("-" * 40)
            self.test_imports()
            self.test_solver_creation()

        if test_type in ["solvers", "all"]:
            print("\n🔧 Testing Solvers")
            print("-" * 40)
            self.test_solver_initialization()
            self.test_solver_information()

        if test_type in ["math", "all"]:
            print("\n🧮 Testing Mathematical Operations")
            print("-" * 40)
            self.test_mathematical_operations()

        if test_type in ["storage", "all"]:
            print("\n💾 Testing Storage Modules")
            print("-" * 40)
            self.test_value_storage()
            self.test_vector_storage()

        if test_type in ["compatibility", "all"]:
            print("\n🔄 Testing Compatibility")
            print("-" * 40)
            self.test_error_handling()
            self.test_memory_management()
            self.test_numerical_stability()

        # Print summary
        print("\n📊 Test Summary")
        print("=" * 60)
        total_tests = len(self.results)
        passed_tests = sum(1 for r in self.results if r.success)
        failed_tests = total_tests - passed_tests

        print(f"Total tests: {total_tests}")
        print(f"Passed: {passed_tests}")
        print(f"Failed: {failed_tests}")

        if self.failed_tests:
            print(f"\n❌ Failed tests: {', '.join(self.failed_tests)}")
            return False
        else:
            print("\n✅ All tests passed!")
            return True


def main():
    """Main function."""
    parser = argparse.ArgumentParser(description="GitHub Actions Test Suite for Komega")
    parser.add_argument(
        "--test-type",
        choices=["basic", "solvers", "math", "storage", "compatibility", "all"],
        default="all",
        help="Type of tests to run",
    )
    parser.add_argument("--verbose", action="store_true", help="Enable verbose output")

    args = parser.parse_args()

    try:
        test_suite = GitHubActionsTestSuite(verbose=args.verbose)
        success = test_suite.run_tests(args.test_type)
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"❌ Test suite failed with error: {e}")
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
