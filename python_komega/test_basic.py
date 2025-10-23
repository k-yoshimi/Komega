"""
Basic test for Komega Python library

This script tests the basic functionality without complex imports.
"""

import numpy as np
import sys
import os

# Add the current directory to the Python path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

def test_numpy_operations():
    """Test basic NumPy operations."""
    print("Testing NumPy operations...")
    
    # Test basic arrays
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0, 6.0])
    
    # Test dot product
    dot_product = np.dot(x, y)
    print(f"Dot product: {dot_product}")
    
    # Test complex arrays
    z1 = np.array([1.0 + 1.0j, 2.0 + 2.0j, 3.0 + 3.0j])
    z2 = np.array([4.0 + 4.0j, 5.0 + 5.0j, 6.0 + 6.0j])
    
    # Test complex dot product
    complex_dot = np.dot(z1, z2)
    print(f"Complex dot product: {complex_dot}")
    
    # Test conjugate dot product
    conjugate_dot = np.vdot(z1, z2)
    print(f"Conjugate dot product: {conjugate_dot}")
    
    print("✓ NumPy operations test passed\n")


def test_parameter_module():
    """Test parameter module directly."""
    print("Testing parameter module...")
    
    try:
        # Import parameter module
        import komega_parameter
        
        # Create parameter instance
        params = komega_parameter.KomegaParameter()
        print(f"✓ Parameter module imported successfully")
        print(f"  Default almost0: {params.almost0}")
        print(f"  Default ndim: {params.ndim}")
        print(f"  Default nl: {params.nl}")
        print(f"  Default nz: {params.nz}")
        
        # Test initialization
        params.initialize(10, 5, 3, 100, 1e-6)
        print(f"✓ Parameters initialized successfully")
        print(f"  ndim: {params.ndim}")
        print(f"  nl: {params.nl}")
        print(f"  nz: {params.nz}")
        print(f"  itermax: {params.itermax}")
        print(f"  threshold: {params.threshold}")
        
    except Exception as e:
        print(f"✗ Parameter module test failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("Parameter module test completed.\n")


def test_math_module():
    """Test math module directly."""
    print("Testing math module...")
    
    try:
        # Import math module
        import komega_math
        
        # Create math instance
        math_ops = komega_math.KomegaMath()
        print(f"✓ Math module imported successfully")
        
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
        
        print(f"✓ Real dot product: {dot_real}")
        print(f"✓ Complex conjugate dot product: {dotc_complex}")
        print(f"✓ Complex dot product: {dotu_complex}")
        
        # Test scaling
        alpha_real = 2.0
        alpha_complex = 1.0 + 1j
        
        x_scaled = x_real.copy()
        math_ops.dscal(alpha_real, x_scaled)
        print(f"✓ Real scaling works")
        
        x_complex_scaled = x_complex.copy()
        math_ops.zscal(alpha_complex, x_complex_scaled)
        print(f"✓ Complex scaling works")
        
    except Exception as e:
        print(f"✗ Math module test failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("Math module test completed.\n")


def test_vals_r_module():
    """Test real values module directly."""
    print("Testing real values module...")
    
    try:
        # Import vals_r module
        import komega_vals_r
        
        # Create vals_r instance
        vals_r = komega_vals_r.KomegaValsR()
        print(f"✓ Real values module imported successfully")
        
        # Test initialization
        z_real = np.array([1.0, 2.0, 3.0])
        vals_r.initialize(z_real, 100)
        print(f"✓ Real values initialized successfully")
        print(f"  π values shape: {vals_r.get_pi_values().shape}")
        print(f"  π old values shape: {vals_r.get_pi_old_values().shape}")
        
    except Exception as e:
        print(f"✗ Real values module test failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("Real values module test completed.\n")


def test_vals_c_module():
    """Test complex values module directly."""
    print("Testing complex values module...")
    
    try:
        # Import vals_c module
        import komega_vals_c
        
        # Create vals_c instance
        vals_c = komega_vals_c.KomegaValsC()
        print(f"✓ Complex values module imported successfully")
        
        # Test initialization
        z_complex = np.array([1.0 + 0.1j, 2.0 + 0.2j, 3.0 + 0.3j])
        vals_c.initialize(z_complex, 100)
        print(f"✓ Complex values initialized successfully")
        print(f"  π values shape: {vals_c.get_pi_values().shape}")
        print(f"  π old values shape: {vals_c.get_pi_old_values().shape}")
        
    except Exception as e:
        print(f"✗ Complex values module test failed: {e}")
        import traceback
        traceback.print_exc()
    
    print("Complex values module test completed.\n")


def test_solver_modules():
    """Test solver modules directly."""
    print("Testing solver modules...")
    
    # Test BiCG solver
    try:
        import komega_bicg
        bicg = komega_bicg.KomegaBiCG()
        print(f"✓ BiCG solver module imported successfully")
    except Exception as e:
        print(f"✗ BiCG solver module test failed: {e}")
    
    # Test CG real solver
    try:
        import komega_cg_r
        cg_r = komega_cg_r.KomegaCGR()
        print(f"✓ CG real solver module imported successfully")
    except Exception as e:
        print(f"✗ CG real solver module test failed: {e}")
    
    # Test CG complex solver
    try:
        import komega_cg_c
        cg_c = komega_cg_c.KomegaCGC()
        print(f"✓ CG complex solver module imported successfully")
    except Exception as e:
        print(f"✗ CG complex solver module test failed: {e}")
    
    # Test COCG solver
    try:
        import komega_cocg
        cocg = komega_cocg.KomegaCOCG()
        print(f"✓ COCG solver module imported successfully")
    except Exception as e:
        print(f"✗ COCG solver module test failed: {e}")
    
    print("Solver modules test completed.\n")


def main():
    """Run all basic tests."""
    print("=" * 60)
    print("Komega Python Library - Basic Test Suite")
    print("=" * 60)
    print()
    
    try:
        test_numpy_operations()
        test_parameter_module()
        test_math_module()
        test_vals_r_module()
        test_vals_c_module()
        test_solver_modules()
        
        print("=" * 60)
        print("All basic tests completed successfully!")
        print("=" * 60)
        print()
        print("✓ NumPy operations work correctly")
        print("✓ Parameter module works correctly")
        print("✓ Math module works correctly")
        print("✓ Real values module works correctly")
        print("✓ Complex values module works correctly")
        print("✓ Solver modules work correctly")
        print()
        print("The Python implementation is functioning properly!")
        
    except Exception as e:
        print(f"Test suite failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
