# Komega C Library

A C port of the Komega Fortran library for solving linear systems using iterative methods.

## Overview

The Komega C library provides a comprehensive set of iterative solvers for linear systems, including:

- **BiCG (Bi-Conjugate Gradient)**: For complex linear systems
- **CG (Conjugate Gradient)**: For real and complex linear systems  
- **COCG (Conjugate Orthogonal Conjugate Gradient)**: For complex linear systems

## Features

- **Modular Design**: Separate modules for parameters, math operations, value storage, and solvers
- **BLAS-like Interface**: Mathematical operations similar to BLAS functions
- **Memory Management**: Automatic memory allocation and cleanup
- **Global Instances**: Compatible with Fortran/Python versions
- **Comprehensive Testing**: Full test suite for verification

## Installation

### Prerequisites

- C compiler (GCC 4.9+ or Clang 3.5+)
- Make
- Math library (libm)

### Building

```bash
# Clone the repository
git clone https://github.com/k-yoshimi/Komega.git
cd Komega/c_komega

# Build the library
make

# Run tests
make run-test

# Install system-wide (optional)
sudo make install
```

### Building Options

```bash
# Build with specific compiler
make CC=gcc

# Build with specific flags
make CFLAGS="-Wall -Wextra -O3"

# Build only static library
make $(LIB_DIR)/libkomega.a

# Build only shared library
make $(LIB_DIR)/libkomega.so
```

## Usage

### Basic Example

```c
#include "komega.h"
#include <stdio.h>

int main() {
    // Create parameter instance
    komega_parameter_t* params = komega_parameter_create();
    komega_parameter_initialize(params, 100, 50, 10, 1000, 1e-6);
    
    // Create math instance
    komega_math_t* math = komega_math_create();
    
    // Test mathematical operations
    double x[5] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double y[5] = {2.0, 3.0, 4.0, 5.0, 6.0};
    double dot_product = komega_math_ddot(x, y, 5);
    printf("Dot product: %f\n", dot_product);
    
    // Cleanup
    komega_math_destroy(math);
    komega_parameter_destroy(params);
    
    return 0;
}
```

### Solver Example

```c
#include "komega.h"
#include <stdio.h>
#include <complex.h>

int main() {
    // Create solver
    komega_solver_t* solver = komega_solver_create(KOMEGA_SOLVER_BICG);
    
    // Initialize with problem parameters
    int ndim = 100, nl = 50, nz = 5;
    int itermax = 1000;
    double threshold = 1e-6;
    
    // Create frequency array
    double complex z[5] = {1.0+0.1*I, 2.0+0.2*I, 3.0+0.3*I, 4.0+0.4*I, 5.0+0.5*I};
    
    // Initialize solver
    komega_status_t status = komega_solver_init(solver, ndim, nl, nz, 
                                               z, itermax, threshold);
    if (status != KOMEGA_SUCCESS) {
        printf("Solver initialization failed: %s\n", 
               komega_get_status_message(status));
        return 1;
    }
    
    // Solve the system (implementation depends on specific solver)
    // ... solver operations ...
    
    // Cleanup
    komega_solver_destroy(solver);
    
    return 0;
}
```

## API Reference

### Core Modules

#### Parameter Module (`komega_parameter.h`)
- `komega_parameter_create()`: Create parameter instance
- `komega_parameter_initialize()`: Initialize with problem parameters
- `komega_parameter_destroy()`: Clean up parameter instance

#### Math Module (`komega_math.h`)
- `komega_math_ddot()`: Real dot product
- `komega_math_zdotc()`: Complex conjugate dot product
- `komega_math_zdotu()`: Complex dot product
- `komega_math_dscal()`: Real scaling
- `komega_math_zscal()`: Complex scaling
- `komega_math_daxpy()`: Real axpy operation
- `komega_math_zaxpy()`: Complex axpy operation

#### Value Storage (`komega_vals_r.h`, `komega_vals_c.h`)
- `komega_vals_r_create()`: Create real value storage
- `komega_vals_c_create()`: Create complex value storage
- `komega_vals_r_initialize()`: Initialize with frequency data
- `komega_vals_c_initialize()`: Initialize with complex frequency data

#### Vector Storage (`komega_vecs_r.h`, `komega_vecs_c.h`)
- `komega_vecs_r_create()`: Create real vector storage
- `komega_vecs_c_create()`: Create complex vector storage
- `komega_vecs_r_initialize()`: Initialize with dimensions
- `komega_vecs_c_initialize()`: Initialize with dimensions

### Solver Interface

#### Solver Creation
```c
komega_solver_t* solver = komega_solver_create(KOMEGA_SOLVER_BICG);
```

#### Solver Initialization
```c
komega_status_t status = komega_solver_init(solver, ndim, nl, nz, 
                                           z_frequencies, itermax, threshold);
```

#### Solver Operations
```c
komega_status_t status = komega_solver_update(solver, v12, v2, x, r_l, status);
```

#### Solver Cleanup
```c
komega_solver_destroy(solver);
```

## Testing

The library includes a comprehensive test suite:

```bash
# Run all tests
make run-test

# Build test executable only
make test

# Run specific test
./build/test_komega
```

### Test Categories

- **Parameter Module**: Parameter creation, initialization, and accessors
- **Math Module**: Mathematical operations (dot products, scaling, axpy, copy)
- **Value Storage**: Real and complex value storage operations
- **Vector Storage**: Real and complex vector storage operations
- **Utility Functions**: Version info, solver names, status messages

## Error Handling

The library uses status codes for error handling:

```c
komega_status_t status = komega_solver_init(solver, ...);
if (status != KOMEGA_SUCCESS) {
    printf("Error: %s\n", komega_get_status_message(status));
    // Handle error
}
```

### Status Codes

- `KOMEGA_SUCCESS`: Operation successful
- `KOMEGA_ERROR_INVALID_PARAM`: Invalid parameter
- `KOMEGA_ERROR_MEMORY_ALLOCATION`: Memory allocation failed
- `KOMEGA_ERROR_SOLVER_NOT_INITIALIZED`: Solver not initialized
- `KOMEGA_ERROR_CONVERGENCE_FAILED`: Convergence failed
- `KOMEGA_ERROR_UNKNOWN_SOLVER_TYPE`: Unknown solver type

## Memory Management

The library provides automatic memory management:

```c
// Create instances
komega_parameter_t* params = komega_parameter_create();
komega_math_t* math = komega_math_create();

// Use instances
// ... operations ...

// Clean up (always call destroy for created instances)
komega_parameter_destroy(params);
komega_math_destroy(math);
```

## Global Instances

For compatibility with Fortran/Python versions, global instances are available:

```c
// Get global instances
komega_parameter_t* params = komega_get_global_parameter();
komega_math_t* math = komega_get_global_math();
komega_vals_r_t* vals_r = komega_get_global_vals_r();
komega_vals_c_t* vals_c = komega_get_global_vals_c();
komega_vecs_r_t* vecs_r = komega_get_global_vecs_r();
komega_vecs_c_t* vecs_c = komega_get_global_vecs_c();

// Clean up global instances
komega_cleanup_global_parameter();
komega_cleanup_global_math();
komega_cleanup_global_vals_r();
komega_cleanup_global_vals_c();
komega_cleanup_global_vecs_r();
komega_cleanup_global_vecs_c();
```

## Compilation Flags

### Recommended Flags

```bash
# Debug build
CFLAGS="-Wall -Wextra -std=c99 -O0 -g -DDEBUG"

# Release build
CFLAGS="-Wall -Wextra -std=c99 -O3 -DNDEBUG"

# With sanitizers
CFLAGS="-Wall -Wextra -std=c99 -O2 -g -fsanitize=address -fsanitize=undefined"
```

### Compiler-Specific Flags

```bash
# GCC
CFLAGS="-Wall -Wextra -std=c99 -O2 -g -fPIC"

# Clang
CFLAGS="-Wall -Wextra -std=c99 -O2 -g -fPIC -Weverything"
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests for new functionality
5. Run the test suite
6. Submit a pull request

## License

This project is licensed under the LGPL v2.1 License - see the LICENSE file for details.

## Acknowledgments

- Original Fortran implementation by Mitsuaki Kawamura
- C port created for verification and testing purposes
- Compatible with Python and Fortran versions

## Version History

- **1.0.0**: Initial C port with basic functionality
- Compatible with Fortran version
- Full test suite included
- GitHub Actions CI/CD support
