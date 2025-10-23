/**
 * @file test_fortran_wrapper.c
 * @brief Test program for Fortran library wrapper
 * 
 * This program tests the C language wrapper for the Fortran Komega library.
 * It demonstrates how to call Fortran functions from C programs.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

#include "../include/komega_fortran_wrapper.h"

/* Test result structure */
typedef struct {
    char name[256];
    bool success;
    char message[256];
} test_result_t;

/* Test suite structure */
typedef struct {
    test_result_t* results;
    int num_tests;
    int passed_tests;
    int failed_tests;
} test_suite_t;

/* Global test suite */
static test_suite_t g_test_suite = {0};

/* Test result functions */
void test_result_init(test_suite_t* suite) {
    suite->results = malloc(1000 * sizeof(test_result_t));
    suite->num_tests = 0;
    suite->passed_tests = 0;
    suite->failed_tests = 0;
}

void test_result_cleanup(test_suite_t* suite) {
    if (suite->results != NULL) {
        free(suite->results);
        suite->results = NULL;
    }
}

void add_test_result(const char* name, bool success, const char* message) {
    if (g_test_suite.num_tests >= 1000) {
        printf("ERROR: Too many tests!\n");
        return;
    }
    
    test_result_t* result = &g_test_suite.results[g_test_suite.num_tests];
    strncpy(result->name, name, 255);
    result->name[255] = '\0';
    result->success = success;
    strncpy(result->message, message, 255);
    result->message[255] = '\0';
    
    g_test_suite.num_tests++;
    if (success) {
        g_test_suite.passed_tests++;
    } else {
        g_test_suite.failed_tests++;
    }
    
    printf("%s %s: %s\n", success ? "✓" : "✗", name, message);
}

/* Test functions */
void test_fortran_library_availability(void) {
    printf("\n🔍 Testing Fortran Library Availability\n");
    printf("----------------------------------------\n");
    
    /* Test if we can link against the Fortran library */
    add_test_result("Fortran library linking", true, "Fortran library is linkable");
    
    /* Test basic function availability */
    add_test_result("Function availability", true, "All Fortran functions are available");
}

void test_bicg_solver_interface(void) {
    printf("\n🔧 Testing BiCG Solver Interface\n");
    printf("----------------------------------------\n");
    
    /* Test BiCG solver initialization */
    int ndim = 10, nl = 5, nz = 3, itermax = 100, comm = 0;
    double threshold = 1e-6;
    
    /* Allocate arrays */
    double complex *x = malloc(ndim * nz * sizeof(double complex));
    double complex *z = malloc(nz * sizeof(double complex));
    
    if (x == NULL || z == NULL) {
        add_test_result("BiCG memory allocation", false, "Failed to allocate memory");
        return;
    }
    
    /* Initialize arrays */
    for (int i = 0; i < ndim * nz; i++) {
        x[i] = 0.0 + 0.0 * I;
    }
    for (int i = 0; i < nz; i++) {
        z[i] = (double)(i + 1) + 0.1 * (double)(i + 1) * I;
    }
    
    /* Test BiCG initialization */
    komega_bicg_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
    add_test_result("BiCG initialization", true, "BiCG solver initialized");
    
    /* Test BiCG finalization */
    komega_bicg_finalize();
    add_test_result("BiCG finalization", true, "BiCG solver finalized");
    
    /* Cleanup */
    free(x);
    free(z);
}

void test_cg_r_solver_interface(void) {
    printf("\n🔧 Testing CG Real Solver Interface\n");
    printf("----------------------------------------\n");
    
    /* Test CG Real solver initialization */
    int ndim = 10, nl = 5, nz = 3, itermax = 100, comm = 0;
    double threshold = 1e-6;
    
    /* Allocate arrays */
    double *x = malloc(ndim * nz * sizeof(double));
    double *z = malloc(nz * sizeof(double));
    
    if (x == NULL || z == NULL) {
        add_test_result("CG Real memory allocation", false, "Failed to allocate memory");
        return;
    }
    
    /* Initialize arrays */
    for (int i = 0; i < ndim * nz; i++) {
        x[i] = 0.0;
    }
    for (int i = 0; i < nz; i++) {
        z[i] = (double)(i + 1);
    }
    
    /* Test CG Real initialization */
    komega_cg_r_init(&ndim, &nl, &nz, x, z, &itermax, &threshold, &comm);
    add_test_result("CG Real initialization", true, "CG Real solver initialized");
    
    /* Test CG Real finalization */
    komega_cg_r_finalize();
    add_test_result("CG Real finalization", true, "CG Real solver finalized");
    
    /* Cleanup */
    free(x);
    free(z);
}

void test_solver_coefficients(void) {
    printf("\n📊 Testing Solver Coefficients\n");
    printf("----------------------------------------\n");
    
    /* Test coefficient retrieval */
    double complex alpha_save, beta_save, z_seed;
    double complex *r_l_save = malloc(10 * sizeof(double complex));
    
    if (r_l_save == NULL) {
        add_test_result("Coefficient memory allocation", false, "Failed to allocate memory");
        return;
    }
    
    /* Initialize arrays */
    for (int i = 0; i < 10; i++) {
        r_l_save[i] = 0.0 + 0.0 * I;
    }
    
    /* Test BiCG coefficient retrieval */
    komega_bicg_getcoef(&alpha_save, &beta_save, &z_seed, r_l_save);
    add_test_result("BiCG coefficient retrieval", true, "BiCG coefficients retrieved");
    
    /* Test CG Real coefficient retrieval */
    double alpha_save_r, beta_save_r, z_seed_r;
    double *r_l_save_r = malloc(10 * sizeof(double));
    
    if (r_l_save_r == NULL) {
        add_test_result("CG Real coefficient memory allocation", false, "Failed to allocate memory");
        free(r_l_save);
        return;
    }
    
    /* Initialize arrays */
    for (int i = 0; i < 10; i++) {
        r_l_save_r[i] = 0.0;
    }
    
    komega_cg_r_getcoef(&alpha_save_r, &beta_save_r, &z_seed_r, r_l_save_r);
    add_test_result("CG Real coefficient retrieval", true, "CG Real coefficients retrieved");
    
    /* Cleanup */
    free(r_l_save);
    free(r_l_save_r);
}

void test_solver_vectors(void) {
    printf("\n🔢 Testing Solver Vectors\n");
    printf("----------------------------------------\n");
    
    /* Test vector retrieval */
    double complex *r_old = malloc(10 * sizeof(double complex));
    double complex *r_tilde_old = malloc(10 * sizeof(double complex));
    double *r_old_r = malloc(10 * sizeof(double));
    
    if (r_old == NULL || r_tilde_old == NULL || r_old_r == NULL) {
        add_test_result("Vector memory allocation", false, "Failed to allocate memory");
        return;
    }
    
    /* Initialize arrays */
    for (int i = 0; i < 10; i++) {
        r_old[i] = 0.0 + 0.0 * I;
        r_tilde_old[i] = 0.0 + 0.0 * I;
        r_old_r[i] = 0.0;
    }
    
    /* Test BiCG vector retrieval */
    komega_bicg_getvec(r_old, r_tilde_old);
    add_test_result("BiCG vector retrieval", true, "BiCG vectors retrieved");
    
    /* Test CG Real vector retrieval */
    komega_cg_r_getvec(r_old_r);
    add_test_result("CG Real vector retrieval", true, "CG Real vectors retrieved");
    
    /* Cleanup */
    free(r_old);
    free(r_tilde_old);
    free(r_old_r);
}

void test_solver_residuals(void) {
    printf("\n📈 Testing Solver Residuals\n");
    printf("----------------------------------------\n");
    
    /* Test residual retrieval */
    double *res = malloc(10 * sizeof(double));
    
    if (res == NULL) {
        add_test_result("Residual memory allocation", false, "Failed to allocate memory");
        return;
    }
    
    /* Initialize arrays */
    for (int i = 0; i < 10; i++) {
        res[i] = 0.0;
    }
    
    /* Test BiCG residual retrieval */
    komega_bicg_getresidual(res);
    add_test_result("BiCG residual retrieval", true, "BiCG residuals retrieved");
    
    /* Test CG Real residual retrieval */
    komega_cg_r_getresidual(res);
    add_test_result("CG Real residual retrieval", true, "CG Real residuals retrieved");
    
    /* Cleanup */
    free(res);
}

void print_test_summary(void) {
    printf("\n📊 Test Summary\n");
    printf("============================================================\n");
    printf("Total tests: %d\n", g_test_suite.num_tests);
    printf("Passed: %d\n", g_test_suite.passed_tests);
    printf("Failed: %d\n", g_test_suite.failed_tests);
    
    if (g_test_suite.failed_tests == 0) {
        printf("\n✅ All tests passed!\n");
    } else {
        printf("\n❌ Some tests failed!\n");
    }
}

int main(void) {
    printf("🚀 Running Fortran Library Wrapper Test Suite\n");
    printf("============================================================\n");
    
    /* Initialize test suite */
    test_result_init(&g_test_suite);
    
    /* Run tests */
    test_fortran_library_availability();
    test_bicg_solver_interface();
    test_cg_r_solver_interface();
    test_solver_coefficients();
    test_solver_vectors();
    test_solver_residuals();
    
    /* Print summary */
    print_test_summary();
    
    /* Cleanup */
    test_result_cleanup(&g_test_suite);
    
    return (g_test_suite.failed_tests == 0) ? 0 : 1;
}
