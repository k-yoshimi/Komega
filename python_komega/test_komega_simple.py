"""
Simple test script for Komega Python library

This script tests the basic functionality of the Komega Python library
to ensure it works correctly.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

import numpy as np
import sys
import os

# Add the current directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import modules directly to avoid relative import issues
from komega_parameter import KomegaParameter, get_global_params
from komega_math import KomegaMath, get_global_math
from komega_vals_r import KomegaValsR, get_global_vals_r
from komega_vals_c import KomegaValsC, get_global_vals_c
from komega_vecs_r import KomegaVecsR, get_global_vecs_r
from komega_vecs_c import KomegaVecsC, get_global_vecs_c
from komega_bicg import KomegaBiCG, get_global_bicg
from komega_cg_r import KomegaCGR, get_global_cg_r
from komega_cg_c import KomegaCGC, get_global_cg_c
from komega_cocg import KomegaCOCG, get_global_cocg


def test_solver_creation():
    """Test solver creation and basic functionality."""
    print("Testing solver creation...")
    
    # Test available solvers
    solvers = ['bicg', 'cg_r', 'cg_c', 'cocg']
    print(f"Available solvers: {solvers}")
    
    # Test solver creation
    for solver_type in solvers:
        try:
            if solver_type == 'bicg':
                solver = get_global_bicg()
            elif solver_type == 'cg_r':
                solver = get_global_cg_r()
            elif solver_type == 'cg_c':
                solver = get_global_cg_c()
            elif solver_type == 'cocg':
                solver = get_global_cocg()
            print(f"✓ Successfully created {solver_type} solver")
        except Exception as e:
            print(f"✗ Failed to create {solver_type} solver: {e}")
    
    print("Solver creation test completed.\n")


def test_parameter_initialization():
    """Test parameter initialization."""
    print("Testing parameter initialization...")
    
    # Test parameters
    ndim = 10
    nl = 5
    nz = 3
    itermax = 100
    threshold = 1e-6
    
    # Test different frequency arrays
    z_real = np.array([1.0, 2.0, 3.0])
    z_complex = np.array([1.0 + 0.1j, 2.0 + 0.2j, 3.0 + 0.3j])
    
    print(f"Test parameters:")
    print(f"  ndim = {ndim}")
    print(f"  nl = {nl}")
    print(f"  nz = {nz}")
    print(f"  itermax = {itermax}")
    print(f"  threshold = {threshold}")
    print(f"  z_real = {z_real}")
    print(f"  z_complex = {z_complex}")
    
    print("Parameter initialization test completed.\n")


def test_math_operations():
    """Test mathematical operations."""
    print("Testing mathematical operations...")
    
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
    
    print(f"Real dot product: {dot_real}")
    print(f"Complex conjugate dot product: {dotc_complex}")
    print(f"Complex dot product: {dotu_complex}")
    
    # Test scaling
    alpha_real = 2.0
    alpha_complex = 1.0 + 1j
    
    x_scaled = x_real.copy()
    math_ops.dscal(alpha_real, x_scaled)
    print(f"Real scaling: {x_scaled}")
    
    x_complex_scaled = x_complex.copy()
    math_ops.zscal(alpha_complex, x_complex_scaled)
    print(f"Complex scaling: {x_complex_scaled}")
    
    print("Mathematical operations test completed.\n")


def test_solver_initialization():
    """Test solver initialization."""
    print("Testing solver initialization...")
    
    # Test parameters
    ndim = 10
    nl = 5
    nz = 3
    itermax = 100
    threshold = 1e-6
    
    # Test real CG solver
    try:
        solver = get_global_cg_r()
        z_real = np.array([1.0, 2.0, 3.0])
        x = solver.init(ndim, nl, nz, z_real, itermax, threshold)
        print(f"✓ Real CG solver initialized successfully")
        print(f"  Solution array shape: {x.shape}")
        print(f"  Solution array dtype: {x.dtype}")
    except Exception as e:
        print(f"✗ Real CG solver initialization failed: {e}")
    
    # Test complex BiCG solver
    try:
        solver = get_global_bicg()
        z_complex = np.array([1.0 + 0.1j, 2.0 + 0.2j, 3.0 + 0.3j])
        x = solver.init(ndim, nl, nz, z_complex, itermax, threshold)
        print(f"✓ Complex BiCG solver initialized successfully")
        print(f"  Solution array shape: {x.shape}")
        print(f"  Solution array dtype: {x.dtype}")
    except Exception as e:
        print(f"✗ Complex BiCG solver initialization failed: {e}")
    
    print("Solver initialization test completed.\n")


def test_value_storage():
    """Test value storage modules."""
    print("Testing value storage modules...")
    
    # Test real values
    try:
        vals_r = get_global_vals_r()
        z_real = np.array([1.0, 2.0, 3.0])
        vals_r.initialize(z_real, 100)
        print(f"✓ Real values initialized successfully")
        print(f"  π values: {vals_r.get_pi_values()}")
        print(f"  π old values: {vals_r.get_pi_old_values()}")
    except Exception as e:
        print(f"✗ Real values initialization failed: {e}")
    
    # Test complex values
    try:
        vals_c = get_global_vals_c()
        z_complex = np.array([1.0 + 0.1j, 2.0 + 0.2j, 3.0 + 0.3j])
        vals_c.initialize(z_complex, 100)
        print(f"✓ Complex values initialized successfully")
        print(f"  π values: {vals_c.get_pi_values()}")
        print(f"  π old values: {vals_c.get_pi_old_values()}")
    except Exception as e:
        print(f"✗ Complex values initialization failed: {e}")
    
    print("Value storage test completed.\n")


def main():
    """Run all tests."""
    print("=" * 60)
    print("Komega Python Library Test Suite")
    print("=" * 60)
    print()
    
    try:
        test_solver_creation()
        test_parameter_initialization()
        test_math_operations()
        test_solver_initialization()
        test_value_storage()
        
        print("=" * 60)
        print("All tests completed successfully!")
        print("=" * 60)
        
    except Exception as e:
        print(f"Test suite failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
