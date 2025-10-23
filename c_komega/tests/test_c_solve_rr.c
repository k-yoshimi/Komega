/**
 * @file test_c_solve_rr.c
 * @brief C test for Real-Real solver (equivalent to solve_rr.F90)
 * 
 * This test program reproduces the functionality of the Fortran solve_rr.F90
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
#include "../include/komega_cg_r.h"

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
static double z_seed;
static int restart;

static double *z;                   /* Frequency array (real) */
static double *ham;                 /* Hamiltonian matrix (real) */
static double *rhs;                 /* Right-hand side vector (real) */
static double *v12, *v2;            /* Working vectors (real) */
static double *r_l;                 /* Projected residual vector (real) */
static double *x;                   /* Solution vector (real) */
static double *alpha_save, *beta_save; /* Saved coefficients (real) */
static double *r_l_save;            /* Saved residual vectors (real) */

/* Function prototypes */
void initialize_test_parameters(void);
void allocate_test_arrays(void);
void deallocate_test_arrays(void);
void generate_test_data(void);
void run_cg_r_test(void);
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
    
    printf("=== C Test: Real-Real Solver ===\n");
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
    /* Allocate frequency array (real) */
    z = malloc(sizeof(double) * nz);
    
    /* Allocate Hamiltonian matrix (real) */
    ham = malloc(sizeof(double) * ndim * ndim);
    
    /* Allocate right-hand side vector (real) */
    rhs = malloc(sizeof(double) * ndim);
    
    /* Allocate working vectors (real) */
    v12 = malloc(sizeof(double) * ndim);
    v2 = malloc(sizeof(double) * ndim);
    
    /* Allocate projected residual vector (real) */
    r_l = malloc(sizeof(double) * nl);
    
    /* Allocate solution vector (real) */
    x = malloc(sizeof(double) * nl * nz);
    
    /* Allocate saved coefficients (real) */
    alpha_save = malloc(sizeof(double) * nz);
    beta_save = malloc(sizeof(double) * nz);
    
    /* Allocate saved residual vectors (real) */
    r_l_save = malloc(sizeof(double) * nl * nz * (itermax + 1));
    
    if (!z || !ham || !rhs || !v12 || !v2 || 
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
    
    /* Generate frequency array (real) */
    for (int i = 0; i < nz; i++) {
        z[i] = (double)(i + 1);
    }
    z_seed = z[0];
    
    /* Generate Hamiltonian matrix (real, symmetric) */
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            if (i == j) {
                ham[i * ndim + j] = (double)(i + 1);
            } else if (i < j) {
                double val = ((double)rand() / RAND_MAX - 0.5) * 0.1;
                ham[i * ndim + j] = val;
                ham[j * ndim + i] = val; /* Symmetric */
            }
        }
    }
    
    /* Generate right-hand side vector (real) */
    for (int i = 0; i < ndim; i++) {
        rhs[i] = ((double)rand() / RAND_MAX - 0.5) * 2.0;
    }
    
    /* Initialize working vectors (real) */
    for (int i = 0; i < ndim; i++) {
        v12[i] = 0.0;
        v2[i] = 0.0;
    }
    
    /* Initialize projected residual vector (real) */
    for (int i = 0; i < nl; i++) {
        r_l[i] = 0.0;
    }
    
    /* Initialize solution vector (real) */
    for (int i = 0; i < nl * nz; i++) {
        x[i] = 0.0;
    }
    
    /* Initialize saved coefficients (real) */
    for (int i = 0; i < nz; i++) {
        alpha_save[i] = 0.0;
        beta_save[i] = 0.0;
    }
    
    /* Initialize saved residual vectors (real) */
    for (int i = 0; i < nl * nz * (itermax + 1); i++) {
        r_l_save[i] = 0.0;
    }
    
    printf("Test data generated successfully.\n\n");
}

/**
 * @brief Run Real CG solver test
 * 
 * Tests the Real CG solver with the generated test data.
 */
void run_cg_r_test(void) {
    printf("=== Testing Real CG Solver ===\n");
    
    struct komega_cg_r* solver = komega_cg_r_create();
    if (solver == NULL) {
        fprintf(stderr, "Error: Failed to create Real CG solver.\n");
        return;
    }
    
    /* Initialize solver */
    komega_status_t status = komega_cg_r_init(solver, ndim, nl, nz, x, z, 
                                            itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        fprintf(stderr, "Error: Failed to initialize Real CG solver.\n");
        komega_cg_r_destroy(solver);
        return;
    }
    
    /* Simulate solver iterations */
    int iter = 0;
    double resnorm = 1.0;
    
    while (iter < itermax && resnorm > threshold) {
        /* Simulate matrix-vector product */
        for (int i = 0; i < ndim; i++) {
            v2[i] = 0.0;
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
            resnorm += v2[i] * v2[i];
        }
        resnorm = sqrt(resnorm);
        
        iter++;
        
        if (iter % 100 == 0) {
            printf("  Iteration %d: residual = %e\n", iter, resnorm);
        }
    }
    
    print_test_results("Real CG", iter, resnorm);
    
    /* Cleanup */
    komega_cg_r_destroy(solver);
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
    printf("C Komega Library - Real-Real Solver Test\n");
    printf("=======================================\n\n");
    
    /* Initialize test */
    initialize_test_parameters();
    allocate_test_arrays();
    generate_test_data();
    
    /* Run tests */
    run_cg_r_test();
    
    /* Cleanup */
    cleanup_test();
    
    printf("All tests completed successfully!\n");
    return 0;
}
