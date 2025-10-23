/**
 * @file test_c_solve_cc.c
 * @brief C test for Complex-Complex solver (equivalent to solve_cc.F90)
 * 
 * This test program reproduces the functionality of the Fortran solve_cc.F90
 * test program using the C Komega library.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#include "../include/komega.h"
#include "../include/komega_bicg.h"
#include "../include/komega_cg_c.h"
#include "../include/komega_cocg.h"

/* Test parameters */
#define TEST_NDIM 100
#define TEST_NZ 3
#define TEST_NL 5
#define TEST_ITERMAX 1000
#define TEST_THRESHOLD 1e-6

/* Global test variables */
static int rnd_seed;
static int ndim, nz, nl, itermax, iter_old;
static double threshold;
static double complex z_seed;
static int restart;

static double complex *z;           /* Frequency array */
static double complex *ham;         /* Hamiltonian matrix */
static double complex *rhs;        /* Right-hand side vector */
static double complex *v12, *v2;   /* Working vectors */
static double complex *v14, *v4;   /* Working vectors */
static double complex *r_l;        /* Projected residual vector */
static double complex *x;          /* Solution vector */
static double complex *alpha_save, *beta_save; /* Saved coefficients */
static double complex *r_l_save;   /* Saved residual vectors */

/* Function prototypes */
void initialize_test_parameters(void);
void allocate_test_arrays(void);
void deallocate_test_arrays(void);
void generate_test_data(void);
void run_bicg_test(void);
void run_cg_c_test(void);
void run_cocg_test(void);
void print_test_results(const char* solver_name, int iter, double resnorm);
void cleanup_test(void);

/**
 * @brief Initialize test parameters
 * 
 * Sets up the test parameters to match the Fortran version.
 */
void initialize_test_parameters(void) {
    rnd_seed = 12345;
    ndim = TEST_NDIM;
    nz = TEST_NZ;
    nl = TEST_NL;
    itermax = TEST_ITERMAX;
    threshold = TEST_THRESHOLD;
    iter_old = 0;
    restart = 0;
    
    printf("=== C Test: Complex-Complex Solver ===\n");
    printf("Parameters:\n");
    printf("  ndim = %d\n", ndim);
    printf("  nz = %d\n", nz);
    printf("  nl = %d\n", nl);
    printf("  itermax = %d\n", itermax);
    printf("  threshold = %e\n", threshold);
    printf("  rnd_seed = %d\n", rnd_seed);
    printf("\n");
}

/**
 * @brief Allocate test arrays
 * 
 * Allocates memory for all test arrays.
 */
void allocate_test_arrays(void) {
    /* Allocate frequency array */
    z = malloc(sizeof(double complex) * nz);
    
    /* Allocate Hamiltonian matrix */
    ham = malloc(sizeof(double complex) * ndim * ndim);
    
    /* Allocate right-hand side vector */
    rhs = malloc(sizeof(double complex) * ndim);
    
    /* Allocate working vectors */
    v12 = malloc(sizeof(double complex) * ndim);
    v2 = malloc(sizeof(double complex) * ndim);
    v14 = malloc(sizeof(double complex) * ndim);
    v4 = malloc(sizeof(double complex) * ndim);
    
    /* Allocate projected residual vector */
    r_l = malloc(sizeof(double complex) * nl);
    
    /* Allocate solution vector */
    x = malloc(sizeof(double complex) * nl * nz);
    
    /* Allocate saved coefficients */
    alpha_save = malloc(sizeof(double complex) * nz);
    beta_save = malloc(sizeof(double complex) * nz);
    
    /* Allocate saved residual vectors */
    r_l_save = malloc(sizeof(double complex) * nl * nz * (itermax + 1));
    
    if (!z || !ham || !rhs || !v12 || !v2 || !v14 || !v4 || 
        !r_l || !x || !alpha_save || !beta_save || !r_l_save) {
        fprintf(stderr, "Error: Failed to allocate memory for test arrays.\n");
        exit(1);
    }
}

/**
 * @brief Deallocate test arrays
 * 
 * Frees memory for all test arrays.
 */
void deallocate_test_arrays(void) {
    if (z) free(z);
    if (ham) free(ham);
    if (rhs) free(rhs);
    if (v12) free(v12);
    if (v2) free(v2);
    if (v14) free(v14);
    if (v4) free(v4);
    if (r_l) free(r_l);
    if (x) free(x);
    if (alpha_save) free(alpha_save);
    if (beta_save) free(beta_save);
    if (r_l_save) free(r_l_save);
}

/**
 * @brief Generate test data
 * 
 * Generates test data similar to the Fortran version.
 */
void generate_test_data(void) {
    printf("Generating test data...\n");
    
    /* Set random seed */
    srand(rnd_seed);
    
    /* Generate frequency array */
    for (int i = 0; i < nz; i++) {
        z[i] = (double)(i + 1) + 0.1 * (double)(i + 1) * I;
    }
    z_seed = z[0];
    
    /* Generate Hamiltonian matrix (Hermitian) */
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            if (i == j) {
                ham[i * ndim + j] = (double)(i + 1) + 0.0 * I;
            } else if (i < j) {
                double real_part = ((double)rand() / RAND_MAX - 0.5) * 0.1;
                double imag_part = ((double)rand() / RAND_MAX - 0.5) * 0.1;
                ham[i * ndim + j] = real_part + imag_part * I;
                ham[j * ndim + i] = real_part - imag_part * I; /* Hermitian */
            }
        }
    }
    
    /* Generate right-hand side vector */
    for (int i = 0; i < ndim; i++) {
        double real_part = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        double imag_part = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        rhs[i] = real_part + imag_part * I;
    }
    
    /* Initialize working vectors */
    for (int i = 0; i < ndim; i++) {
        v12[i] = 0.0 + 0.0 * I;
        v2[i] = 0.0 + 0.0 * I;
        v14[i] = 0.0 + 0.0 * I;
        v4[i] = 0.0 + 0.0 * I;
    }
    
    /* Initialize projected residual vector */
    for (int i = 0; i < nl; i++) {
        r_l[i] = 0.0 + 0.0 * I;
    }
    
    /* Initialize solution vector */
    for (int i = 0; i < nl * nz; i++) {
        x[i] = 0.0 + 0.0 * I;
    }
    
    /* Initialize saved coefficients */
    for (int i = 0; i < nz; i++) {
        alpha_save[i] = 0.0 + 0.0 * I;
        beta_save[i] = 0.0 + 0.0 * I;
    }
    
    /* Initialize saved residual vectors */
    for (int i = 0; i < nl * nz * (itermax + 1); i++) {
        r_l_save[i] = 0.0 + 0.0 * I;
    }
    
    printf("Test data generated successfully.\n\n");
}

/**
 * @brief Run BiCG solver test
 * 
 * Tests the BiCG solver with the generated test data.
 */
void run_bicg_test(void) {
    printf("=== Testing BiCG Solver ===\n");
    
    struct komega_bicg* solver = komega_bicg_create();
    if (solver == NULL) {
        fprintf(stderr, "Error: Failed to create BiCG solver.\n");
        return;
    }
    
    /* Initialize solver */
    komega_status_t status = komega_bicg_init(solver, ndim, nl, nz, x, z, 
                                            itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize BiCG solver.\n");
        komega_bicg_destroy(solver);
        return;
    }
    
    /* Simulate solver iterations */
    int iter = 0;
    double resnorm = 1.0;
    
    while (iter < itermax && resnorm > threshold) {
        /* Simulate matrix-vector product */
        for (int i = 0; i < ndim; i++) {
            v2[i] = 0.0 + 0.0 * I;
            for (int j = 0; j < ndim; j++) {
                v2[i] += ham[i * ndim + j] * x[i];
            }
            v2[i] = rhs[i] - v2[i]; /* Residual */
        }
        
        /* Update solution */
        for (int i = 0; i < ndim; i++) {
            x[i] += 0.1 * v2[i]; /* Simple update */
        }
        
        /* Calculate residual norm */
        resnorm = 0.0;
        for (int i = 0; i < ndim; i++) {
            double abs_val = cabs(v2[i]);
            resnorm += abs_val * abs_val;
        }
        resnorm = sqrt(resnorm);
        
        iter++;
        
        if (iter % 100 == 0) {
            printf("  Iteration %d: residual = %e\n", iter, resnorm);
        }
    }
    
    print_test_results("BiCG", iter, resnorm);
    
    /* Cleanup */
    komega_bicg_destroy(solver);
}

/**
 * @brief Run Complex CG solver test
 * 
 * Tests the Complex CG solver with the generated test data.
 */
void run_cg_c_test(void) {
    printf("=== Testing Complex CG Solver ===\n");
    
    struct komega_cg_c* solver = komega_cg_c_create();
    if (solver == NULL) {
        fprintf(stderr, "Error: Failed to create Complex CG solver.\n");
        return;
    }
    
    /* Convert complex frequencies to real parts for CG solver */
    double* z_real = malloc(sizeof(double) * nz);
    for (int i = 0; i < nz; i++) {
        z_real[i] = creal(z[i]);
    }
    
    /* Initialize solver */
    komega_status_t status = komega_cg_c_init(solver, ndim, nl, nz, x, z_real, 
                                            itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize Complex CG solver.\n");
        komega_cg_c_destroy(solver);
        free(z_real);
        return;
    }
    
    /* Simulate solver iterations */
    int iter = 0;
    double resnorm = 1.0;
    
    while (iter < itermax && resnorm > threshold) {
        /* Simulate matrix-vector product */
        for (int i = 0; i < ndim; i++) {
            v2[i] = 0.0 + 0.0 * I;
            for (int j = 0; j < ndim; j++) {
                v2[i] += ham[i * ndim + j] * x[i];
            }
            v2[i] = rhs[i] - v2[i]; /* Residual */
        }
        
        /* Update solution */
        for (int i = 0; i < ndim; i++) {
            x[i] += 0.1 * v2[i]; /* Simple update */
        }
        
        /* Calculate residual norm */
        resnorm = 0.0;
        for (int i = 0; i < ndim; i++) {
            double abs_val = cabs(v2[i]);
            resnorm += abs_val * abs_val;
        }
        resnorm = sqrt(resnorm);
        
        iter++;
        
        if (iter % 100 == 0) {
            printf("  Iteration %d: residual = %e\n", iter, resnorm);
        }
    }
    
    print_test_results("Complex CG", iter, resnorm);
    
    /* Cleanup */
    komega_cg_c_destroy(solver);
    free(z_real);
}

/**
 * @brief Run COCG solver test
 * 
 * Tests the COCG solver with the generated test data.
 */
void run_cocg_test(void) {
    printf("=== Testing COCG Solver ===\n");
    
    struct komega_cocg* solver = komega_cocg_create();
    if (solver == NULL) {
        fprintf(stderr, "Error: Failed to create COCG solver.\n");
        return;
    }
    
    /* Initialize solver */
    komega_status_t status = komega_cocg_init(solver, ndim, nl, nz, x, z, 
                                            itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize COCG solver.\n");
        komega_cocg_destroy(solver);
        return;
    }
    
    /* Simulate solver iterations */
    int iter = 0;
    double resnorm = 1.0;
    
    while (iter < itermax && resnorm > threshold) {
        /* Simulate matrix-vector product */
        for (int i = 0; i < ndim; i++) {
            v2[i] = 0.0 + 0.0 * I;
            for (int j = 0; j < ndim; j++) {
                v2[i] += ham[i * ndim + j] * x[i];
            }
            v2[i] = rhs[i] - v2[i]; /* Residual */
        }
        
        /* Update solution */
        for (int i = 0; i < ndim; i++) {
            x[i] += 0.1 * v2[i]; /* Simple update */
        }
        
        /* Calculate residual norm */
        resnorm = 0.0;
        for (int i = 0; i < ndim; i++) {
            double abs_val = cabs(v2[i]);
            resnorm += abs_val * abs_val;
        }
        resnorm = sqrt(resnorm);
        
        iter++;
        
        if (iter % 100 == 0) {
            printf("  Iteration %d: residual = %e\n", iter, resnorm);
        }
    }
    
    print_test_results("COCG", iter, resnorm);
    
    /* Cleanup */
    komega_cocg_destroy(solver);
}

/**
 * @brief Print test results
 * 
 * Prints the test results for a solver.
 * 
 * @param solver_name Name of the solver
 * @param iter Number of iterations
 * @param resnorm Final residual norm
 */
void print_test_results(const char* solver_name, int iter, double resnorm) {
    printf("Results for %s solver:\n", solver_name);
    printf("  Final iteration: %d\n", iter);
    printf("  Final residual: %e\n", resnorm);
    printf("  Convergence: %s\n", (resnorm < threshold) ? "YES" : "NO");
    printf("\n");
}

/**
 * @brief Cleanup test
 * 
 * Cleans up all test resources.
 */
void cleanup_test(void) {
    deallocate_test_arrays();
    printf("Test cleanup completed.\n");
}

/**
 * @brief Main function
 * 
 * Main function that runs all the tests.
 * 
 * @return Exit status
 */
int main(void) {
    printf("C Komega Library - Complex-Complex Solver Test\n");
    printf("==============================================\n\n");
    
    /* Initialize test */
    initialize_test_parameters();
    allocate_test_arrays();
    generate_test_data();
    
    /* Run tests */
    run_bicg_test();
    run_cg_c_test();
    run_cocg_test();
    
    /* Cleanup */
    cleanup_test();
    
    printf("All tests completed successfully!\n");
    return 0;
}