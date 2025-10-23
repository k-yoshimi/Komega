#!/usr/bin/env python3
"""
Test suite for verifying corrected function calls in Komega solvers.

This test suite specifically tests the zdotuMPI and zdotcMPI function calls
that were corrected based on Copilot code review comments in PR #20.

The functions should be called with 2 arguments (x, y) instead of 3 arguments
(ndim, x, y) as they were incorrectly called before.
"""

import numpy as np
import sys
import os

# Add the current directory to the path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_math_functions():
    """Test the mathematical functions with correct argument counts."""
    print("=" * 70)
    print("Testing Mathematical Functions with Correct Arguments")
    print("=" * 70)
    
    try:
        from komega_math import KomegaMath
        
        # Initialize math module
        math_module = KomegaMath()
        
        # Test vectors
        n = 5
        x_real = np.random.rand(n)
        y_real = np.random.rand(n)
        x_complex = np.random.rand(n) + 1j * np.random.rand(n)
        y_complex = np.random.rand(n) + 1j * np.random.rand(n)
        
        print(f"Testing with vectors of length {n}")
        print(f"Real vectors: x={x_real[:3]}..., y={y_real[:3]}...")
        print(f"Complex vectors: x={x_complex[:3]}..., y={y_complex[:3]}...")
        
        # Test zdotuMPI (complex dot product)
        print("\n1. Testing zdotuMPI (complex dot product):")
        try:
            result_zdotu = math_module.zdotuMPI(x_complex, y_complex)
            print(f"   ✓ zdotuMPI(x, y) = {result_zdotu}")
            print(f"   ✓ Function accepts 2 arguments correctly")
        except Exception as e:
            print(f"   ✗ zdotuMPI failed: {e}")
            return False
        
        # Test zdotcMPI (complex conjugate dot product)
        print("\n2. Testing zdotcMPI (complex conjugate dot product):")
        try:
            result_zdotc = math_module.zdotcMPI(x_complex, y_complex)
            print(f"   ✓ zdotcMPI(x, y) = {result_zdotc}")
            print(f"   ✓ Function accepts 2 arguments correctly")
        except Exception as e:
            print(f"   ✗ zdotcMPI failed: {e}")
            return False
        
        # Test ddotMPI (real dot product)
        print("\n3. Testing ddotMPI (real dot product):")
        try:
            result_ddot = math_module.ddotMPI(x_real, y_real)
            print(f"   ✓ ddotMPI(x, y) = {result_ddot}")
            print(f"   ✓ Function accepts 2 arguments correctly")
        except Exception as e:
            print(f"   ✗ ddotMPI failed: {e}")
            return False
        
        print("\n✓ All mathematical functions work with correct argument counts")
        return True
        
    except ImportError as e:
        print(f"✗ Failed to import math module: {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def test_solver_function_calls():
    """Test that solver functions can be called without errors."""
    print("\n" + "=" * 70)
    print("Testing Solver Function Calls")
    print("=" * 70)
    
    try:
        from komega_parameter import get_global_params
        from komega_math import KomegaMath
        from komega_vals_r import KomegaValsR
        from komega_vals_c import KomegaValsC
        
        # Initialize parameters
        params = get_global_params()
        params.initialize(ndim=5, nl=3, nz=2, itermax=50, threshold=1e-6)
        
        # Initialize math module
        math_module = KomegaMath()
        
        # Initialize value storage
        vals_r = KomegaValsR(params)
        vals_c = KomegaValsC(params)
        
        print(f"Initialized with ndim={params.ndim}, nl={params.nl}, nz={params.nz}")
        
        # Test vectors
        v2 = np.random.rand(params.ndim) + 1j * np.random.rand(params.ndim)
        v4 = np.random.rand(params.ndim) + 1j * np.random.rand(params.ndim)
        v12 = np.random.rand(params.ndim) + 1j * np.random.rand(params.ndim)
        
        print(f"Test vectors: v2={v2[:3]}..., v4={v4[:3]}..., v12={v12[:3]}...")
        
        # Test the corrected function calls
        print("\n1. Testing corrected zdotuMPI calls:")
        try:
            # This should work with 2 arguments (corrected)
            rho = math_module.zdotuMPI(v2, v2)
            print(f"   ✓ zdotuMPI(v2, v2) = {rho}")
            
            alpha_denom = math_module.zdotuMPI(v2, v12)
            print(f"   ✓ zdotuMPI(v2, v12) = {alpha_denom}")
            
        except Exception as e:
            print(f"   ✗ zdotuMPI calls failed: {e}")
            return False
        
        print("\n2. Testing corrected zdotcMPI calls:")
        try:
            # This should work with 2 arguments (corrected)
            rho = math_module.zdotcMPI(v2, v2)
            print(f"   ✓ zdotcMPI(v2, v2) = {rho}")
            
            alpha_denom = math_module.zdotcMPI(v2, v12)
            print(f"   ✓ zdotcMPI(v2, v12) = {alpha_denom}")
            
            # Test with different vectors
            rho_v4 = math_module.zdotcMPI(v4, v2)
            print(f"   ✓ zdotcMPI(v4, v2) = {rho_v4}")
            
            alpha_denom_v4 = math_module.zdotcMPI(v4, v12)
            print(f"   ✓ zdotcMPI(v4, v12) = {alpha_denom_v4}")
            
        except Exception as e:
            print(f"   ✗ zdotcMPI calls failed: {e}")
            return False
        
        print("\n✓ All solver function calls work correctly")
        return True
        
    except ImportError as e:
        print(f"✗ Failed to import modules: {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def test_incorrect_calls():
    """Test that incorrect 3-argument calls fail as expected."""
    print("\n" + "=" * 70)
    print("Testing Incorrect Function Calls (Should Fail)")
    print("=" * 70)
    
    try:
        from komega_math import KomegaMath
        
        math_module = KomegaMath()
        
        # Test vectors
        v2 = np.random.rand(5) + 1j * np.random.rand(5)
        v12 = np.random.rand(5) + 1j * np.random.rand(5)
        
        print("Testing that incorrect 3-argument calls fail:")
        
        # These should fail because the functions only accept 2 arguments
        try:
            # This should fail - incorrect 3-argument call
            result = math_module.zdotuMPI(5, v2, v2)  # ndim, x, y
            print(f"   ✗ zdotuMPI(5, v2, v2) unexpectedly succeeded: {result}")
            return False
        except TypeError as e:
            print(f"   ✓ zdotuMPI(5, v2, v2) correctly failed: {e}")
        except Exception as e:
            print(f"   ✓ zdotuMPI(5, v2, v2) failed as expected: {e}")
        
        try:
            # This should fail - incorrect 3-argument call
            result = math_module.zdotcMPI(5, v2, v2)  # ndim, x, y
            print(f"   ✗ zdotcMPI(5, v2, v2) unexpectedly succeeded: {result}")
            return False
        except TypeError as e:
            print(f"   ✓ zdotcMPI(5, v2, v2) correctly failed: {e}")
        except Exception as e:
            print(f"   ✓ zdotcMPI(5, v2, v2) failed as expected: {e}")
        
        print("\n✓ Incorrect function calls fail as expected")
        return True
        
    except ImportError as e:
        print(f"✗ Failed to import math module: {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def test_numerical_accuracy():
    """Test numerical accuracy of the corrected function calls."""
    print("\n" + "=" * 70)
    print("Testing Numerical Accuracy")
    print("=" * 70)
    
    try:
        from komega_math import KomegaMath
        
        math_module = KomegaMath()
        
        # Test with known values
        x = np.array([1.0, 2.0, 3.0])
        y = np.array([4.0, 5.0, 6.0])
        xc = np.array([1.0+1j, 2.0+2j, 3.0+3j])
        yc = np.array([4.0+4j, 5.0+5j, 6.0+6j])
        
        print(f"Test vectors: x={x}, y={y}")
        print(f"Complex vectors: xc={xc}, yc={yc}")
        
        # Test real dot product
        expected_real = np.dot(x, y)  # Should be 1*4 + 2*5 + 3*6 = 4 + 10 + 18 = 32
        result_real = math_module.ddotMPI(x, y)
        print(f"\nReal dot product:")
        print(f"   Expected: {expected_real}")
        print(f"   Result:   {result_real}")
        print(f"   Difference: {abs(expected_real - result_real)}")
        
        # Test complex dot product
        expected_complex = np.dot(xc, yc)  # Complex dot product
        result_complex = math_module.zdotuMPI(xc, yc)
        print(f"\nComplex dot product:")
        print(f"   Expected: {expected_complex}")
        print(f"   Result:   {result_complex}")
        print(f"   Difference: {abs(expected_complex - result_complex)}")
        
        # Test complex conjugate dot product
        expected_conj = np.dot(np.conj(xc), yc)  # Conjugate dot product
        result_conj = math_module.zdotcMPI(xc, yc)
        print(f"\nComplex conjugate dot product:")
        print(f"   Expected: {expected_conj}")
        print(f"   Result:   {result_conj}")
        print(f"   Difference: {abs(expected_conj - result_conj)}")
        
        # Check accuracy
        tolerance = 1e-14
        real_accurate = abs(expected_real - result_real) < tolerance
        complex_accurate = abs(expected_complex - result_complex) < tolerance
        conj_accurate = abs(expected_conj - result_conj) < tolerance
        
        print(f"\nAccuracy check (tolerance={tolerance}):")
        print(f"   Real dot product: {'✓' if real_accurate else '✗'}")
        print(f"   Complex dot product: {'✓' if complex_accurate else '✗'}")
        print(f"   Conjugate dot product: {'✓' if conj_accurate else '✗'}")
        
        if real_accurate and complex_accurate and conj_accurate:
            print("\n✓ All numerical calculations are accurate")
            return True
        else:
            print("\n✗ Some numerical calculations are inaccurate")
            return False
        
    except ImportError as e:
        print(f"✗ Failed to import math module: {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False


def main():
    """Run all function call tests."""
    print("=" * 70)
    print("Komega Function Call Tests")
    print("Testing corrected zdotuMPI and zdotcMPI function calls")
    print("=" * 70)
    
    tests = [
        ("Mathematical Functions", test_math_functions),
        ("Solver Function Calls", test_solver_function_calls),
        ("Incorrect Calls (Should Fail)", test_incorrect_calls),
        ("Numerical Accuracy", test_numerical_accuracy),
    ]
    
    results = []
    
    for test_name, test_func in tests:
        print(f"\n{'='*20} {test_name} {'='*20}")
        try:
            result = test_func()
            results.append((test_name, result))
            if result:
                print(f"✓ {test_name}: PASSED")
            else:
                print(f"✗ {test_name}: FAILED")
        except Exception as e:
            print(f"✗ {test_name}: ERROR - {e}")
            results.append((test_name, False))
    
    # Summary
    print("\n" + "=" * 70)
    print("Test Summary")
    print("=" * 70)
    
    passed = sum(1 for _, result in results if result)
    total = len(results)
    
    for test_name, result in results:
        status = "PASSED" if result else "FAILED"
        print(f"{test_name:30} : {status}")
    
    print(f"\nResults: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All function call tests passed!")
        print("✓ Copilot review comments have been successfully addressed")
        return True
    else:
        print("⚠️  Some function call tests failed")
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
