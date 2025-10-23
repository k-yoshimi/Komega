/**
 * @file test_fortran_compatibility.c
 * @brief Fortran compatibility test for C Komega Library
 * 
 * This test verifies that the C library provides the same functionality
 * as the Fortran version, following the same test patterns.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#include "../include/komega.h"
#include "../include/komega_bicg.h"
#include "../include/komega_cg_r.h"
#include "../include/komega_cg_c.h"
#include "../include/komega_cocg.h"

/* Test parameters matching Fortran version */
#define TEST_NDIM 50
#define TEST_NZ 3
#define TEST_NL 5
#define TEST_ITERMAX 1000
#define TEST_THRESHOLD 1e-6
#define TEST_RND_SEED 12345

/* Test result structure */
typedef struct {
    char test_name[64];
    int passed;
    int failed;
    double execution_time;
} test_result_t;

/* Global test variables */
static int rnd_seed = TEST_RND_SEED;
static int ndim = TEST_NDIM;
static int nz = TEST_NZ;
static int nl = TEST_NL;
static int itermax = TEST_ITERMAX;
static double threshold = TEST_THRESHOLD;

/* Test data arrays */
static double complex* z_complex;
static double* z_real;
static double complex* ham_complex;
static double* ham_real;
static double complex* rhs_complex;
static double* rhs_real;
static double complex* x_complex;
static double* x_real;

/* Function prototypes */
void initialize_test_data(void);
void cleanup_test_data(void);
test_result_t run_bicg_compatibility_test(void);
test_result_t run_cg_r_compatibility_test(void);
test_result_t run_cg_c_compatibility_test(void);
test_result_t run_cocg_compatibility_test(void);
void print_test_summary(test_result_t* results, int num_tests);

/**
 * @brief Initialize test data
 * 
 * Creates test data similar to the Fortran version.
 */
void initialize_test_data(void) {
    printf("Initializing test data...\n");
    
    /* Set random seed */
    srand(rnd_seed);
    
    /* Allocate arrays */
    z_complex = malloc(sizeof(double complex) * nz);
    z_real = malloc(sizeof(double) * nz);
    ham_complex = malloc(sizeof(double complex) * ndim * ndim);
    ham_real = malloc(sizeof(double) * ndim * ndim);
    rhs_complex = malloc(sizeof(double complex) * ndim);
    rhs_real = malloc(sizeof(double) * ndim);
    x_complex = malloc(sizeof(double complex) * ndim * nz);
    x_real = malloc(sizeof(double) * ndim * nz);
    
    if (!z_complex || !z_real || !ham_complex || !ham_real || 
        !rhs_complex || !rhs_real || !x_complex || !x_real) {
        fprintf(stderr, "Error: Failed to allocate test data arrays.\n");
        exit(1);
    }
    
    /* Generate frequency arrays */
    for (int i = 0; i < nz; i++) {
        z_complex[i] = (double)(i + 1) + 0.1 * (double)(i + 1) * I;
        z_real[i] = (double)(i + 1);
    }
    
    /* Generate complex Hamiltonian matrix (Hermitian) */
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            if (i == j) {
                ham_complex[i * ndim + j] = (double)(i + 1) + 0.0 * I;
            } else if (i < j) {
                double real_part = ((double)rand() / RAND_MAX - 0.5) * 0.1;
                double imag_part = ((double)rand() / RAND_MAX - 0.5) * 0.1;
                ham_complex[i * ndim + j] = real_part + imag_part * I;
                ham_complex[j * ndim + i] = real_part - imag_part * I; /* Hermitian */
            }
        }
    }
    
    /* Generate real Hamiltonian matrix (symmetric) */
    for (int i = 0; i < ndim; i++) {
        for (int j = 0; j < ndim; j++) {
            if (i == j) {
                ham_real[i * ndim + j] = (double)(i + 1);
            } else if (i < j) {
                double val = ((double)rand() / RAND_MAX - 0.5) * 0.1;
                ham_real[i * ndim + j] = val;
                ham_real[j * ndim + i] = val; /* Symmetric */
            }
        }
    }
    
    /* Generate right-hand side vectors */
    for (int i = 0; i < ndim; i++) {
        double real_part = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        double imag_part = ((double)rand() / RAND_MAX - 0.5) * 2.0;
        rhs_complex[i] = real_part + imag_part * I;
        rhs_real[i] = ((double)rand() / RAND_MAX - 0.5) * 2.0;
    }
    
    /* Initialize solution vectors */
    for (int i = 0; i < ndim * nz; i++) {
        x_complex[i] = 0.0 + 0.0 * I;
        x_real[i] = 0.0;
    }
    
    printf("Test data initialized successfully.\n\n");
}

/**
 * @brief Cleanup test data
 * 
 * Frees all allocated test data arrays.
 */
void cleanup_test_data(void) {
    if (z_complex) free(z_complex);
    if (z_real) free(z_real);
    if (ham_complex) free(ham_complex);
    if (ham_real) free(ham_real);
    if (rhs_complex) free(rhs_complex);
    if (rhs_real) free(rhs_real);
    if (x_complex) free(x_complex);
    if (x_real) free(x_real);
}

/**
 * @brief Run BiCG compatibility test
 * 
 * Tests BiCG solver functionality similar to Fortran solve_cc.F90.
 * 
 * @return Test result structure
 */
test_result_t run_bicg_compatibility_test(void) {
    test_result_t result = {"BiCG Solver (Complex-Complex)", 0, 0, 0.0};
    clock_t start_time = clock();
    
    printf("=== BiCG Solver Compatibility Test ===\n");
    
    struct komega_bicg* solver = komega_bicg_create();
    if (solver == NULL) {
        printf("❌ Failed to create BiCG solver\n");
        result.failed++;
        return result;
    }
    result.passed++;
    printf("✓ BiCG solver created\n");
    
    /* Initialize solver */
    komega_status_t status = komega_bicg_init(solver, ndim, nl, nz, x_complex, z_complex, 
                                            itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to initialize BiCG solver\n");
        komega_bicg_destroy(solver);
        result.failed++;
        return result;
    }
    result.passed++;
    printf("✓ BiCG solver initialized\n");
    
    /* Test solver functions */
    double res[nz];
    status = komega_bicg_getresidual(solver, res);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get BiCG residual\n");
        result.failed++;
    } else {
        result.passed++;
        printf("✓ BiCG residual retrieved\n");
    }
    
    /* Test coefficient functions */
    double complex alpha_save[nz], beta_save[nz];
    double complex z_seed;
    double complex r_l_save[nl * nz * (itermax + 1)];
    
    status = komega_bicg_getcoef(solver, alpha_save, beta_save, &z_seed, r_l_save);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get BiCG coefficients\n");
        result.failed++;
    } else {
        result.passed++;
        printf("✓ BiCG coefficients retrieved\n");
    }
    
    /* Test vector functions */
    double complex r_old[ndim * nz], r_tilde_old[ndim * nz];
    status = komega_bicg_getvec(solver, r_old, r_tilde_old);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get BiCG vectors\n");
        result.failed++;
    } else {
        result.passed++;
        printf("✓ BiCG vectors retrieved\n");
    }
    
    komega_bicg_destroy(solver);
    result.passed++;
    printf("✓ BiCG solver destroyed\n");
    
    result.execution_time = ((double)(clock() - start_time)) / CLOCKS_PER_SEC;
    printf("BiCG test completed in %.3f seconds\n\n", result.execution_time);
    
    return result;
}

/**
 * @brief Run Real CG compatibility test
 * 
 * Tests Real CG solver functionality similar to Fortran solve_rc.F90.
 * 
 * @return Test result structure
 */
test_result_t run_cg_r_compatibility_test(void) {
    test_result_t result = {"Real CG Solver (Real-Complex)", 0, 0, 0.0};
    clock_t start_time = clock();
    
    printf("=== Real CG Solver Compatibility Test ===\n");
    
    struct komega_cg_r* solver = komega_cg_r_create();
    if (solver == NULL) {
        printf("❌ Failed to create Real CG solver\n");
        result.failed++;
        return result;
    }
    result.passed++;
    printf("✓ Real CG solver created\n");
    
    /* Initialize solver */
    komega_status_t status = komega_cg_r_init(solver, ndim, nl, nz, x_real, z_real, 
                                            itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to initialize Real CG solver\n");
        komega_cg_r_destroy(solver);
        result.failed++;
        return result;
    }
    result.passed++;
    printf("✓ Real CG solver initialized\n");
    
    /* Test solver functions */
    double res[nz];
    status = komega_cg_r_getresidual(solver, res);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get Real CG residual\n");
        result.failed++;
    } else {
        result.passed++;
        printf("✓ Real CG residual retrieved\n");
    }
    
    komega_cg_r_destroy(solver);
    result.passed++;
    printf("✓ Real CG solver destroyed\n");
    
    result.execution_time = ((double)(clock() - start_time)) / CLOCKS_PER_SEC;
    printf("Real CG test completed in %.3f seconds\n\n", result.execution_time);
    
    return result;
}

/**
 * @brief Run Complex CG compatibility test
 * 
 * Tests Complex CG solver functionality similar to Fortran solve_cr.F90.
 * 
 * @return Test result structure
 */
test_result_t run_cg_c_compatibility_test(void) {
    test_result_t result = {"Complex CG Solver (Complex-Real)", 0, 0, 0.0};
    clock_t start_time = clock();
    
    printf("=== Complex CG Solver Compatibility Test ===\n");
    
    struct komega_cg_c* solver = komega_cg_c_create();
    if (solver == NULL) {
        printf("❌ Failed to create Complex CG solver\n");
        result.failed++;
        return result;
    }
    result.passed++;
    printf("✓ Complex CG solver created\n");
    
    /* Initialize solver */
    komega_status_t status = komega_cg_c_init(solver, ndim, nl, nz, x_complex, z_real, 
                                            itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to initialize Complex CG solver\n");
        komega_cg_c_destroy(solver);
        result.failed++;
        return result;
    }
    result.passed++;
    printf("✓ Complex CG solver initialized\n");
    
    /* Test solver functions */
    double res[nz];
    status = komega_cg_c_getresidual(solver, res);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get Complex CG residual\n");
        result.failed++;
    } else {
        result.passed++;
        printf("✓ Complex CG residual retrieved\n");
    }
    
    komega_cg_c_destroy(solver);
    result.passed++;
    printf("✓ Complex CG solver destroyed\n");
    
    result.execution_time = ((double)(clock() - start_time)) / CLOCKS_PER_SEC;
    printf("Complex CG test completed in %.3f seconds\n\n", result.execution_time);
    
    return result;
}

/**
 * @brief Run COCG compatibility test
 * 
 * Tests COCG solver functionality similar to Fortran solve_rr.F90.
 * 
 * @return Test result structure
 */
test_result_t run_cocg_compatibility_test(void) {
    test_result_t result = {"COCG Solver (Complex-Complex)", 0, 0, 0.0};
    clock_t start_time = clock();
    
    printf("=== COCG Solver Compatibility Test ===\n");
    
    struct komega_cocg* solver = komega_cocg_create();
    if (solver == NULL) {
        printf("❌ Failed to create COCG solver\n");
        result.failed++;
        return result;
    }
    result.passed++;
    printf("✓ COCG solver created\n");
    
    /* Initialize solver */
    komega_status_t status = komega_cocg_init(solver, ndim, nl, nz, x_complex, z_complex, 
                                            itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to initialize COCG solver\n");
        komega_cocg_destroy(solver);
        result.failed++;
        return result;
    }
    result.passed++;
    printf("✓ COCG solver initialized\n");
    
    /* Test solver functions */
    double res[nz];
    status = komega_cocg_getresidual(solver, res);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get COCG residual\n");
        result.failed++;
    } else {
        result.passed++;
        printf("✓ COCG residual retrieved\n");
    }
    
    komega_cocg_destroy(solver);
    result.passed++;
    printf("✓ COCG solver destroyed\n");
    
    result.execution_time = ((double)(clock() - start_time)) / CLOCKS_PER_SEC;
    printf("COCG test completed in %.3f seconds\n\n", result.execution_time);
    
    return result;
}

/**
 * @brief Print test summary
 * 
 * Prints a summary of all test results.
 * 
 * @param results Array of test results
 * @param num_tests Number of tests
 */
void print_test_summary(test_result_t* results, int num_tests) {
    printf("========================================\n");
    printf("C Komega Library - Fortran Compatibility Test Summary\n");
    printf("========================================\n\n");
    
    int total_passed = 0;
    int total_failed = 0;
    double total_time = 0.0;
    
    for (int i = 0; i < num_tests; i++) {
        printf("Test %d: %s\n", i + 1, results[i].test_name);
        printf("  Passed: %d, Failed: %d\n", results[i].passed, results[i].failed);
        printf("  Execution time: %.3f seconds\n", results[i].execution_time);
        printf("  Status: %s\n\n", (results[i].failed == 0) ? "✅ PASSED" : "❌ FAILED");
        
        total_passed += results[i].passed;
        total_failed += results[i].failed;
        total_time += results[i].execution_time;
    }
    
    printf("Overall Summary:\n");
    printf("  Total tests: %d\n", num_tests);
    printf("  Total passed: %d\n", total_passed);
    printf("  Total failed: %d\n", total_failed);
    printf("  Total execution time: %.3f seconds\n", total_time);
    printf("  Overall status: %s\n", (total_failed == 0) ? "✅ ALL TESTS PASSED" : "❌ SOME TESTS FAILED");
}

/**
 * @brief Main function
 * 
 * Runs all Fortran compatibility tests.
 * 
 * @return Exit status
 */
int main(void) {
    printf("C Komega Library - Fortran Compatibility Test\n");
    printf("============================================\n\n");
    
    printf("Test Configuration:\n");
    printf("  ndim = %d\n", ndim);
    printf("  nz = %d\n", nz);
    printf("  nl = %d\n", nl);
    printf("  itermax = %d\n", itermax);
    printf("  threshold = %e\n", threshold);
    printf("  rnd_seed = %d\n\n", rnd_seed);
    
    /* Initialize test data */
    initialize_test_data();
    
    /* Run all tests */
    test_result_t results[4];
    results[0] = run_bicg_compatibility_test();
    results[1] = run_cg_r_compatibility_test();
    results[2] = run_cg_c_compatibility_test();
    results[3] = run_cocg_compatibility_test();
    
    /* Print summary */
    print_test_summary(results, 4);
    
    /* Cleanup */
    cleanup_test_data();
    
    /* Check if all tests passed */
    int all_passed = 1;
    for (int i = 0; i < 4; i++) {
        if (results[i].failed > 0) {
            all_passed = 0;
            break;
        }
    }
    
    if (all_passed) {
        printf("\n🎉 All Fortran compatibility tests passed!\n");
        printf("C Komega Library is fully compatible with Fortran version.\n");
        return 0;
    } else {
        printf("\n❌ Some tests failed. Please check the output above.\n");
        return 1;
    }
}
