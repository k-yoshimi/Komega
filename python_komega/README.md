# Komega Python Library

A Python port of the Komega Fortran library for solving linear systems in materials science.

## Overview

This library provides Python implementations of the Komega Fortran library, including:

- **BiCG (Bi-Conjugate Gradient)** solver for complex linear systems
- **CG (Conjugate Gradient)** solver for real and complex linear systems  
- **COCG (Conjugate Orthogonal Conjugate Gradient)** solver for complex linear systems
- **QMR (Quasi-Minimal Residual)** solver for complex linear systems

## Features

- Multiple solver algorithms for different problem types
- Support for multiple frequency points (shifted equations)
- MPI support for distributed computing
- Restart capability for long-running calculations
- Unified interface for easy solver switching
- Comprehensive test suite

## Installation

The Python library is located in the `python_komega` directory. To use it:

```python
import sys
sys.path.append('/path/to/python_komega')
from komega import create_solver, get_available_solvers
```

## Quick Start

### Basic Usage

```python
import numpy as np
from komega import create_solver

# Create a BiCG solver
solver = create_solver('bicg')

# Initialize the solver
ndim = 100      # Dimension of Hamiltonian matrix
nl = 50         # Dimension of projection space
nz = 10         # Number of frequency points
z = np.array([1.0 + 0.1j, 2.0 + 0.2j, 3.0 + 0.3j, ...])  # Frequencies
itermax = 1000  # Maximum iterations
threshold = 1e-6  # Convergence threshold

x = solver.init(ndim, nl, nz, z, itermax, threshold)

# Update the solver (in a loop)
status = [0, 0, 0]
solver.update(v12, v2, v14, v4, x, r_l, status)

# Check convergence
if status[0] < 0:
    print(f"Converged in {abs(status[0])} iterations")
else:
    print(f"Continuing iteration {status[0]}")

# Get residual norms
residuals = solver.get_residual()
print(f"Residual norms: {residuals}")

# Clean up
solver.finalize()
```

### Available Solvers

```python
from komega import get_available_solvers, get_solver_info

# Get available solvers
solvers = get_available_solvers()
print(f"Available solvers: {solvers}")

# Get solver information
info = get_solver_info('bicg')
print(f"BiCG solver: {info['name']} - {info['description']}")
```

### Solver Types

| Solver | Type | Description | Data Type |
|--------|------|-------------|-----------|
| `bicg` | BiCG | Bi-Conjugate Gradient | Complex |
| `cg_r` | CG Real | Conjugate Gradient | Real |
| `cg_c` | CG Complex | Conjugate Gradient | Complex |
| `cocg` | COCG | Conjugate Orthogonal CG | Complex |

## API Reference

### Main Interface

#### `create_solver(solver_type: str) -> KomegaSolver`
Create a new solver instance.

**Parameters:**
- `solver_type`: Type of solver ('bicg', 'cg_r', 'cg_c', 'cocg')

**Returns:**
- `KomegaSolver`: New solver instance

#### `get_available_solvers() -> List[str]`
Get list of available solver types.

**Returns:**
- `List[str]`: List of available solver types

#### `get_solver_info(solver_type: str) -> dict`
Get information about a specific solver.

**Parameters:**
- `solver_type`: Type of solver

**Returns:**
- `dict`: Solver information

### Solver Methods

#### `solver.init(ndim, nl, nz, z, itermax, threshold, comm=None) -> np.ndarray`
Initialize the solver.

**Parameters:**
- `ndim`: Dimension of the Hamiltonian matrix
- `nl`: Dimension of the projection space
- `nz`: Number of frequency points
- `z`: Frequency array
- `itermax`: Maximum number of iterations
- `threshold`: Convergence threshold
- `comm`: MPI communicator (optional)

**Returns:**
- `np.ndarray`: Initialized solution array

#### `solver.update(*args, **kwargs) -> None`
Update the solver iteration.

**Parameters:**
- Vary depending on solver type

#### `solver.get_coefficients() -> tuple`
Get saved coefficients for restart.

**Returns:**
- `tuple`: (alpha_save, beta_save, z_seed, r_l_save)

#### `solver.get_vectors() -> np.ndarray or tuple`
Get old residual vectors.

**Returns:**
- `np.ndarray` or `tuple`: Old residual vector(s)

#### `solver.get_residual() -> np.ndarray`
Get residual norms for all frequencies.

**Returns:**
- `np.ndarray`: Residual norms

#### `solver.finalize() -> None`
Finalize the solver and clean up memory.

## Testing

The library includes comprehensive test suites to verify functionality and compatibility with the original Fortran implementation.

### Quick Test

Run the basic test suite:

```python
python test_komega.py
```

### Complete Test Suite

Run all tests including Fortran compatibility tests:

```bash
python run_tests.py
```

### Test Types

- **Basic Tests**: Core functionality and API tests
- **Fortran Compatibility**: Tests that reproduce Fortran test cases
- **Detailed Compatibility**: Comprehensive comparison with Fortran results

### Individual Test Suites

```bash
# Basic functionality tests
python run_tests.py --test-type basic

# Fortran compatibility tests  
python run_tests.py --test-type fortran

# Detailed Fortran compatibility tests
python run_tests.py --test-type detailed

# All tests (default)
python run_tests.py --test-type all
```

### Test Output

The test suites provide detailed output including:
- Solver initialization and configuration
- Iteration progress and convergence status
- Residual norms and solution vectors
- Comparison with expected Fortran results

## Examples

### Example 1: Basic BiCG Solver

```python
import numpy as np
from komega import create_solver

# Create solver
solver = create_solver('bicg')

# Problem parameters
ndim = 100
nl = 50
nz = 5
z = np.array([1.0 + 0.1j, 2.0 + 0.2j, 3.0 + 0.3j, 4.0 + 0.4j, 5.0 + 0.5j])
itermax = 1000
threshold = 1e-6

# Initialize
x = solver.init(ndim, nl, nz, z, itermax, threshold)

# Solve (pseudo-code - actual implementation would depend on your problem)
for iteration in range(itermax):
    # Compute residual and other vectors
    # ... (problem-specific code)
    
    # Update solver
    status = [0, 0, 0]
    solver.update(v12, v2, v14, v4, x, r_l, status)
    
    # Check convergence
    if status[0] < 0:
        print(f"Converged in {abs(status[0])} iterations")
        break
    elif status[0] == itermax:
        print("Maximum iterations reached")
        break

# Get results
residuals = solver.get_residual()
print(f"Final residuals: {residuals}")

# Clean up
solver.finalize()
```

### Example 2: Solver Comparison

```python
import numpy as np
from komega import create_solver, get_available_solvers

# Problem parameters
ndim = 100
nl = 50
nz = 3
z = np.array([1.0 + 0.1j, 2.0 + 0.2j, 3.0 + 0.3j])
itermax = 1000
threshold = 1e-6

# Test different solvers
solvers_to_test = ['bicg', 'cg_c', 'cocg']

for solver_type in solvers_to_test:
    print(f"Testing {solver_type} solver...")
    
    try:
        solver = create_solver(solver_type)
        x = solver.init(ndim, nl, nz, z, itermax, threshold)
        print(f"  ✓ {solver_type} initialized successfully")
        solver.finalize()
    except Exception as e:
        print(f"  ✗ {solver_type} failed: {e}")
```

## Notes

- This is a Python port of the original Fortran library
- The API is designed to be compatible with the Fortran version
- MPI support is available but requires proper MPI setup
- The library is designed for materials science applications
- All solvers support restart capability for long-running calculations

## License

This library is licensed under the GNU Lesser General Public License v3.0.

## Contact

For questions or issues, please contact the original author:
- **Author**: Kazuyoshi Yoshimi
- **Email**: k-yoshimi@issp.u-tokyo.ac.jp
