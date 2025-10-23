# Test Suite for Komega Library

This directory contains comprehensive tests for the Komega linear equation solver library. The tests verify the functionality and performance of various iterative solvers for different matrix-vector combinations.

## Test Programs

### 1. solve_rr (Real-Real Solver)
- **Purpose**: Tests real matrix-vector linear equation solver
- **Algorithm**: Conjugate Gradient (CG) method
- **Input**: Real frequency values from `real_freq.in`
- **Test Coverage**:
  - Real symmetric positive definite matrix generation
  - CG convergence behavior
  - Residual vector computation and verification
  - Restart functionality testing

### 2. solve_rc (Real-Complex Solver)
- **Purpose**: Tests real matrix with complex vector linear equation solver
- **Algorithm**: Conjugate Orthogonal CG (COCG) method
- **Input**: Complex frequency values from `complex_freq.in`
- **Test Coverage**:
  - Real matrix with complex vector solving
  - COCG convergence behavior
  - Complex residual vector computation and verification
  - Orthogonality verification

### 3. solve_cr (Complex-Real Solver)
- **Purpose**: Tests complex matrix with real vector linear equation solver
- **Algorithm**: Complex CG (CG_C) method
- **Input**: Real frequency values from `real_freq.in`
- **Test Coverage**:
  - Complex matrix generation (Hermitian)
  - CG_C convergence behavior
  - Real vector with complex matrix operations

### 4. solve_cc (Complex-Complex Solver)
- **Purpose**: Tests complex matrix with complex vector linear equation solver
- **Algorithm**: Bi-Conjugate Gradient (BiCG) method
- **Input**: Complex frequency values from `complex_freq.in`
- **Test Coverage**:
  - Complex matrix generation (Hermitian)
  - BiCG convergence behavior
  - Complex residual vector computation and verification
  - Biorthogonality verification

## Supporting Modules

### mathlib.F90
- Interface definitions for BLAS/LAPACK routines
- Provides standardized mathematical operations for all test programs

### diagonalize.F90
- Diagonalization test module (incomplete)
- Intended for eigenvalue problem testing

### make_ham.F90
- Utility for generating random positive definite matrices
- Used for creating test matrices with known properties

## Input Files

### real_freq.in
- Contains real frequency values for testing
- Format: Real numbers (e.g., -2.0, -1.0, 0.0, 1.0, 2.0)

### complex_freq.in
- Contains complex frequency values for testing
- Format: Complex numbers (e.g., (-2.0, 1.0), (-1.0, 1.0), (0.0, 1.0))

## Test Features

### Convergence Testing
- Verifies convergence behavior of each algorithm
- Tests with different matrix sizes and frequency values
- Monitors iteration counts and convergence criteria

### Accuracy Testing
- Computes residual vectors to verify solution accuracy
- Compares against expected results
- Tests numerical stability

### Restart Functionality
- Tests computation interruption and resumption
- Verifies state preservation across restarts
- Ensures consistency in restart scenarios

### Orthogonality Testing
- Verifies orthogonality of iteration vectors
- Tests biorthogonality for BiCG method
- Ensures numerical stability of iterative processes

### Memory Management
- Tests proper allocation and deallocation of dynamic arrays
- Verifies memory efficiency
- Tests for memory leaks

## Building and Running Tests

### Prerequisites
- Fortran 90/95 compiler (gfortran recommended)
- BLAS and LAPACK libraries
- Make build system

### Build Commands
```bash
# Build all test programs
make

# Build specific test
make solve_rr.x
make solve_rc.x
make solve_cr.x
make solve_cc.x
```

### Running Tests
```bash
# Run real-real test
./solve_rr.x < real_freq.in

# Run real-complex test
./solve_rc.x < complex_freq.in

# Run complex-real test
./solve_cr.x < real_freq.in

# Run complex-complex test
./solve_cc.x < complex_freq.in
```

## Test Parameters

Each test program accepts the following parameters via namelist input:

- `ndim`: Dimension of the vector space
- `nl`: Number of left vectors
- `nz`: Number of frequencies
- `itermax`: Maximum number of iterations
- `threshold`: Convergence threshold
- `rnd_seed`: Random number generator seed
- `restart`: Restart flag for continuation runs

## Expected Output

Each test program provides:
- Input parameter summary
- Linear system generation details
- Iteration progress and convergence status
- Final solution vectors
- Residual vector verification
- Restart file generation (if applicable)

## Troubleshooting

### Common Issues
1. **Convergence failures**: Check matrix conditioning and threshold settings
2. **Memory errors**: Verify BLAS/LAPACK library linkage
3. **Compilation errors**: Ensure proper Fortran compiler and library paths

### Debug Information
- Each test provides detailed debug output during execution
- Status codes indicate convergence behavior
- Residual norms are printed for accuracy verification

## Contributing

When adding new tests:
1. Follow the existing naming convention (solve_xx format)
2. Include comprehensive test coverage
3. Add appropriate input files
4. Update this README with test descriptions
5. Ensure proper error handling and cleanup
