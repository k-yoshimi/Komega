<a name= "english">

[日本語](README_ja.md) / English 

<img src="doc/figs/komega.png" width="300">

## Build Status
[![CI](https://github.com/k-yoshimi/Komega/workflows/CI/badge.svg)](https://github.com/k-yoshimi/Komega/actions)
[![Extended Test Matrix](https://github.com/k-yoshimi/Komega/workflows/Extended%20Test%20Matrix/badge.svg)](https://github.com/k-yoshimi/Komega/actions)
[![Quick Test](https://github.com/k-yoshimi/Komega/workflows/Quick%20Test/badge.svg)](https://github.com/k-yoshimi/Komega/actions)

## Project Information
[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg)](https://www.gnu.org/licenses/lgpl-3.0)
[![Language: Fortran](https://img.shields.io/badge/Language-Fortran-orange.svg)](https://fortran-lang.org/)
[![Platform: Linux, macOS](https://img.shields.io/badge/Platform-Linux%20%7C%20macOS-lightgrey.svg)](https://github.com/k-yoshimi/Komega)
[![Version: 2.0.0](https://img.shields.io/badge/Version-2.0.0-green.svg)](https://github.com/k-yoshimi/Komega/releases)

# What is Kω? 

This package provides the solver library based on Shifted-Krylov subspace method and the software to calculate dynamical Green function by inputting the Hamiltonian and the vector of the excited state.

# Download

 * Latest release

   https://github.com/issp-center-dev/Komega/releases/download/v2.0.0/komega-2.0.0.tar.gz
 * Release note

   https://github.com/issp-center-dev/Komega/releases/tag/v2.0.0
 * Other version
 
   https://github.com/issp-center-dev/Komega/releases
   
# Prerequisite

 * Fortran compiler (gfortran 10, 11, or 12 recommended)
 * BLAS library (reference BLAS or OpenBLAS)
 * LAPACK library (Used in the sample program)
 * MPI library (Optional)
 * Autotools (autoconf, automake, libtool) for building from source

# Docments

 * Manual for the Library
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/library/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/ja/_build/latex/komega.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/library/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/library/en/_build/latex/komega.pdf))
 * Manual for the sample program
   * Japanese ([HTML](https://issp-center-dev.github.io/Komega/software/ja/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/ja/_build/latex/shiftk.pdf))
   * English ([HTML](https://issp-center-dev.github.io/Komega/software/en/_build/html/index.html)/[PDF](https://issp-center-dev.github.io/Komega/software/en/_build/latex/shiftk.pdf))

# License

GNU LESSER GENERAL PUBLIC LICENSE Version 3

The development of Kω is supported by “Project for advancement of software usability in materials science” of ISSP, The University of Tokyo.

(＊) We hope that you cite the following paper and repository when you publish the results using Kω.

Papaer: : [“Kω — Open-source library for the shifted Krylov subspace method of the form (zI−H)x=b”, Takeo Hoshi, Mitsuaki Kawamura, Kazuyoshi Yoshimi, Yuichi Motoyama, Takahiro Misawa, Youhei Yamaji, Synge Todo, Naoki Kawashima, Tomohiro Sogabe, Computer Physics Communications, Volume 258, January 2021, 107536.](https://www.sciencedirect.com/science/article/pii/S0010465520302551)
Repository: https://github.com/issp-center-dev/Komega



# Directory Tree

 * `app/`: The directory for the software
   * `src/`: The source directory for the software
   * `sample/`: The sample directory
     * `Shiftk.nb`: The mathmatica note book to test the software (for developers)
     * `denovo/`: Sample files not to input either Hamiltonian and the right-hand side initial vector
     * `from_file/`: Sample files to input both Hamiltonian and the right-hand side initial vector
 * `doc/`: The documents directory (only japanese, english version will be provided from ver.1.0)
   * `index.html` : Top page for documents
   * `library/`: Directory for the document of the librariy
   * `software/`: Directory for the document of the sample program
 * `configure`: Configuration script to build
 * `src/`: The source directory for the libraries
 * `test/`: The test directory for the libraries

# Install

## Quick Start

The simplest procedure is as follows:

```bash
$ ./configure --prefix=install_dir
$ make
$ make install
```

where `install_dir` should be replaced with the full path of the directory where
the library will be stored.

## Generated Files

The following objects are generated in the directory specified by `install_dir`:
 * In `install_dir/lib/`: Static and shared libraries (`libkomega.a`, `libkomega.so`)
 * In `install_dir/include/`: Header file for C/C++ (`komega.h`)
 * `install_dir/bin/ShiftK.out`: Sample program

## Build Options

```bash
# Configure with OpenBLAS
$ ./configure --prefix=install_dir --with-blas=openblas

# Configure with reference BLAS (default)
$ ./configure --prefix=install_dir
```

For more details, please see the manual.

# Test of the software

## Sample Applications

### Denovo Sample (No input files required)
```bash
$ cd app/sample/denovo/
$ install_dir/bin/ShiftK.out namelist.def lBiCG
```

### From File Sample (With input files)
```bash
$ cd app/sample/from_file/
$ install_dir/bin/ShiftK.out namelist.def lBiCG
```

## Expected Output

When the software works well, the following files will be generated:
 * `output/dynamicalG.dat`: Dynamical Green's function data
 * `output/ResVec.dat0`: Residual vector data
 * `output/TriDiagComp.dat`: Tridiagonal components data

The details of the file format of `namelist.def` is written in the manual.

# Usage of libraries

## How to call each routines in the program

### For fortran/C/C++

See the manual.

## How to link the libraries

```
$ gfortran myprog.f90 -L install_dir/lib -lkomega -lblas -I install_dir/include
$ gcc myprog.c -L install_dir/lib -lkomega -lblas -I install_dir/include
```
etc.

Add the `install_dir/lib` directry to the environment
variable `LD_LIBRARY_PATH` to execute the file with dynamic link.

# Test for libraries (Optional)

## Running Library Tests

Move to the `test/` directory and run the comprehensive test suite:

```bash
$ cd test/
$ make
$ ./run_tests.sh
```

## Individual Solver Tests

You can test individual solvers:

```bash
# Complex-Complex solver
$ ./solve_cc.x < complex_freq.in

# Real-Complex solver  
$ ./solve_rc.x < complex_freq.in

# Complex-Real solver
$ ./solve_cr.x < real_freq.in

# Real-Real solver
$ ./solve_rr.x < real_freq.in

# QMR solvers
$ ./solve_rc_qmr1.x < real_freq.in
$ ./solve_rc_qmr2.x < real_freq.in
```

## Test Parameters

The parameters in `krylov.in` (you can modify the file name freely) are as follows:
 * `ndim`: The dimension of the pseudo Hamiltonian
 * `nl`: This parameter is used to test the projection.
   The vectors are calculated from the target vector up to `nl(<=ndim)`th vector.
 * `nz`: The number of the frequencies to calculate.
 * `itermax`: The maximum number of iterations.
 * `threshold`: The threshold to judge the convergence.
 * `rnd_seed`: The seed of random number to generate pseudo Hamiltonian.
 * Write the value of each frequencies line by line after this namelist.

# Continuous Integration

This project uses GitHub Actions for automated testing across multiple environments:

## Tested Environments
 * **Ubuntu 22.04**: gfortran-10 (reference BLAS), gfortran-11 (OpenBLAS)
 * **Ubuntu Latest**: gfortran-12 (OpenBLAS)

## Test Coverage
 * **Sample Applications**: Denovo and from_file samples
 * **Library Tests**: All solver algorithms (BiCG, COCG, CG, QMR)
 * **Compiler Compatibility**: Multiple gfortran versions
 * **BLAS Libraries**: Both reference BLAS and OpenBLAS
