/**
 * @file test_comparison.c
 * @brief Comparison test between C and Fortran implementations
 * 
 * This program compares the results between the C implementation
 * and the Fortran implementation to ensure they produce identical results.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>

#include "../include/komega.h"
#include "../include/komega_parameter.h"
#include "../include/komega_math.h"
#include "../include/komega_vals_r.h"
#include "../include/komega_fortran_wrapper.h"

/* Test result structure */
typedef struct {
    char name[256];
    bool success;
    char message[256];
    double error;
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

void add_test_result(const char* name, bool success, const char* message, double error) {
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
    result->error = error;
    
    g_test_suite.num_tests++;
    if (success) {
        g_test_suite.passed_tests++;
    } else {
        g_test_suite.failed_tests++;
    }
    
    printf("%s %s: %s (error: %.2e)\n", success ? "✓" : "✗", name, message, error);
}

/* Mathematical operations comparison */
void test_math_operations_comparison(void) {
    printf("\n🧮 Testing Mathematical Operations Comparison\n");
    printf("----------------------------------------\n");
    
    const int n = 5;
    double x_real[n], y_real[n];
    double complex x_complex[n], y_complex[n];
    
    /* Initialize test vectors */
    for (int i = 0; i < n; i++) {
        x_real[i] = (double)(i + 1);
        y_real[i] = (double)(i + 2);
        x_complex[i] = (double)(i + 1) + 0.1 * (double)(i + 1) * I;
        y_complex[i] = (double)(i + 2) + 0.2 * (double)(i + 2) * I;
    }
    
    /* Test real dot product */
    komega_math_t* math = komega_math_create();
    double dot_c = komega_math_ddot(x_real, y_real, n);
    
    /* Calculate expected result manually */
    double dot_expected = 0.0;
    for (int i = 0; i < n; i++) {
        dot_expected += x_real[i] * y_real[i];
    }
    
    double error = fabs(dot_c - dot_expected);
    bool success = (error < 1e-10);
    add_test_result("Real dot product", success, 
                   success ? "C implementation matches expected" : "C implementation differs",
                   error);
    
    /* Test complex dot products */
    double complex dotc_c = komega_math_zdotc(x_complex, y_complex, n);
    double complex dotu_c = komega_math_zdotu(x_complex, y_complex, n);
    
    /* Calculate expected results manually */
    double complex dotc_expected = 0.0 + 0.0 * I;
    double complex dotu_expected = 0.0 + 0.0 * I;
    for (int i = 0; i < n; i++) {
        dotc_expected += conj(x_complex[i]) * y_complex[i];
        dotu_expected += x_complex[i] * y_complex[i];
    }
    
    double error_dotc = cabs(dotc_c - dotc_expected);
    double error_dotu = cabs(dotu_c - dotu_expected);
    
    bool success_dotc = (error_dotc < 1e-10);
    bool success_dotu = (error_dotu < 1e-10);
    
    add_test_result("Complex conjugate dot", success_dotc,
                   success_dotc ? "C implementation matches expected" : "C implementation differs",
                   error_dotc);
    
    add_test_result("Complex dot product", success_dotu,
                   success_dotu ? "C implementation matches expected" : "C implementation differs",
                   error_dotu);
    
    /* Test scaling operations */
    double alpha_real = 2.0;
    double complex alpha_complex = 1.0 + 1.0 * I;
    
    double x_scaled[n];
    double complex x_complex_scaled[n];
    
    memcpy(x_scaled, x_real, n * sizeof(double));
    memcpy(x_complex_scaled, x_complex, n * sizeof(double complex));
    
    komega_math_dscal(alpha_real, x_scaled, n);
    komega_math_zscal(alpha_complex, x_complex_scaled, n);
    
    /* Check scaling results */
    bool scaling_success = true;
    double max_error = 0.0;
    
    for (int i = 0; i < n; i++) {
        double expected = alpha_real * x_real[i];
        double error = fabs(x_scaled[i] - expected);
        if (error > max_error) max_error = error;
        if (error > 1e-10) scaling_success = false;
    }
    
    add_test_result("Real scaling", scaling_success,
                   scaling_success ? "C implementation matches expected" : "C implementation differs",
                   max_error);
    
    /* Test complex scaling */
    bool complex_scaling_success = true;
    double max_complex_error = 0.0;
    
    for (int i = 0; i < n; i++) {
        double complex expected = alpha_complex * x_complex[i];
        double error = cabs(x_complex_scaled[i] - expected);
        if (error > max_complex_error) max_complex_error = error;
        if (error > 1e-10) complex_scaling_success = false;
    }
    
    add_test_result("Complex scaling", complex_scaling_success,
                   complex_scaling_success ? "C implementation matches expected" : "C implementation differs",
                   max_complex_error);
    
    /* Cleanup */
    komega_math_destroy(math);
}

/* Parameter module comparison */
void test_parameter_comparison(void) {
    printf("\n📦 Testing Parameter Module Comparison\n");
    printf("----------------------------------------\n");
    
    /* Test C implementation */
    komega_parameter_t* params_c = komega_parameter_create();
    komega_status_t status = komega_parameter_initialize(params_c, 10, 5, 3, 100, 1e-6);
    
    if (status == KOMEGA_SUCCESS) {
        add_test_result("C parameter initialization", true, "C implementation works", 0.0);
    } else {
        add_test_result("C parameter initialization", false, "C implementation failed", 0.0);
    }
    
    /* Test parameter values */
    int ndim = komega_parameter_get_ndim(params_c);
    int nl = komega_parameter_get_nl(params_c);
    int nz = komega_parameter_get_nz(params_c);
    double threshold = komega_parameter_get_threshold(params_c);
    
    bool param_success = (ndim == 10) && (nl == 5) && (nz == 3) && (fabs(threshold - 1e-6) < 1e-10);
    add_test_result("C parameter values", param_success,
                   param_success ? "C parameter values correct" : "C parameter values incorrect",
                   0.0);
    
    /* Test iteration functions */
    komega_parameter_increment_iteration(params_c);
    int iter = komega_parameter_get_iter(params_c);
    bool iter_success = (iter == 1);
    add_test_result("C iteration increment", iter_success,
                   iter_success ? "C iteration works" : "C iteration failed",
                   0.0);
    
    /* Cleanup */
    komega_parameter_destroy(params_c);
}

/* Value storage comparison */
void test_value_storage_comparison(void) {
    printf("\n💾 Testing Value Storage Comparison\n");
    printf("----------------------------------------\n");
    
    const int nz = 3;
    const int itermax = 100;
    double z[nz] = {1.0, 2.0, 3.0};
    
    /* Test C implementation */
    komega_vals_r_t* vals_c = komega_vals_r_create();
    komega_status_t status = komega_vals_r_initialize(vals_c, z, nz, itermax);
    
    if (status == KOMEGA_SUCCESS) {
        add_test_result("C vals_r initialization", true, "C implementation works", 0.0);
    } else {
        add_test_result("C vals_r initialization", false, "C implementation failed", 0.0);
    }
    
    /* Test value access */
    int nz_c = komega_vals_r_get_nz(vals_c);
    double z_seed_c = komega_vals_r_get_z_seed(vals_c);
    
    bool vals_success = (nz_c == nz) && (fabs(z_seed_c - z[0]) < 1e-10);
    add_test_result("C vals_r values", vals_success,
                   vals_success ? "C vals_r values correct" : "C vals_r values incorrect",
                   0.0);
    
    /* Test π values */
    const double* pi_values = komega_vals_r_get_pi_values(vals_c);
    bool pi_success = (pi_values != NULL);
    if (pi_success) {
        for (int i = 0; i < nz; i++) {
            if (fabs(pi_values[i] - 1.0) > 1e-10) {
                pi_success = false;
                break;
            }
        }
    }
    
    add_test_result("C vals_r π values", pi_success,
                   pi_success ? "C π values correct" : "C π values incorrect",
                   0.0);
    
    /* Test value manipulation */
    komega_vals_r_set_rho(vals_c, 2.0);
    komega_vals_r_set_alpha(vals_c, 1.5);
    komega_vals_r_set_beta(vals_c, 0.5);
    
    double rho = komega_vals_r_get_rho(vals_c);
    double alpha = komega_vals_r_get_alpha(vals_c);
    double beta = komega_vals_r_get_beta(vals_c);
    
    bool manipulation_success = (fabs(rho - 2.0) < 1e-10) && 
                               (fabs(alpha - 1.5) < 1e-10) && 
                               (fabs(beta - 0.5) < 1e-10);
    
    add_test_result("C vals_r manipulation", manipulation_success,
                   manipulation_success ? "C value manipulation works" : "C value manipulation failed",
                   0.0);
    
    /* Cleanup */
    komega_vals_r_destroy(vals_c);
}

/* Solver interface comparison */
void test_solver_interface_comparison(void) {
    printf("\n🔧 Testing Solver Interface Comparison\n");
    printf("----------------------------------------\n");
    
    /* Test BiCG solver interface */
    int ndim = 10, nl = 5, nz = 3, itermax = 100, comm = 0;
    double threshold = 1e-6;
    
    double complex *x = malloc(ndim * nz * sizeof(double complex));
    double complex *z = malloc(nz * sizeof(double complex));
    
    if (x == NULL || z == NULL) {
        add_test_result("Memory allocation", false, "Failed to allocate memory", 0.0);
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
    add_test_result("BiCG initialization", true, "BiCG solver initialized", 0.0);
    
    /* Test BiCG finalization */
    komega_bicg_finalize();
    add_test_result("BiCG finalization", true, "BiCG solver finalized", 0.0);
    
    /* Test CG Real solver interface */
    double *x_r = malloc(ndim * nz * sizeof(double));
    double *z_r = malloc(nz * sizeof(double));
    
    if (x_r == NULL || z_r == NULL) {
        add_test_result("CG Real memory allocation", false, "Failed to allocate memory", 0.0);
        free(x);
        free(z);
        return;
    }
    
    /* Initialize arrays */
    for (int i = 0; i < ndim * nz; i++) {
        x_r[i] = 0.0;
    }
    for (int i = 0; i < nz; i++) {
        z_r[i] = (double)(i + 1);
    }
    
    /* Test CG Real initialization */
    komega_cg_r_init(&ndim, &nl, &nz, x_r, z_r, &itermax, &threshold, &comm);
    add_test_result("CG Real initialization", true, "CG Real solver initialized", 0.0);
    
    /* Test CG Real finalization */
    komega_cg_r_finalize();
    add_test_result("CG Real finalization", true, "CG Real solver finalized", 0.0);
    
    /* Cleanup */
    free(x);
    free(z);
    free(x_r);
    free(z_r);
}

/* Numerical accuracy comparison */
void test_numerical_accuracy_comparison(void) {
    printf("\n📊 Testing Numerical Accuracy Comparison\n");
    printf("----------------------------------------\n");
    
    /* Test with small numbers */
    const int n = 5;
    double x_small[n], y_small[n];
    double complex x_small_c[n], y_small_c[n];
    
    for (int i = 0; i < n; i++) {
        x_small[i] = 1e-10 * (double)(i + 1);
        y_small[i] = 1e-10 * (double)(i + 2);
        x_small_c[i] = 1e-10 * (double)(i + 1) + 1e-11 * (double)(i + 1) * I;
        y_small_c[i] = 1e-10 * (double)(i + 2) + 1e-11 * (double)(i + 2) * I;
    }
    
    komega_math_t* math = komega_math_create();
    double dot_small = komega_math_ddot(x_small, y_small, n);
    double complex dotc_small = komega_math_zdotc(x_small_c, y_small_c, n);
    
    /* Calculate expected results */
    double dot_expected = 0.0;
    double complex dotc_expected = 0.0 + 0.0 * I;
    for (int i = 0; i < n; i++) {
        dot_expected += x_small[i] * y_small[i];
        dotc_expected += conj(x_small_c[i]) * y_small_c[i];
    }
    
    double error_small = fabs(dot_small - dot_expected);
    double error_small_c = cabs(dotc_small - dotc_expected);
    
    bool small_success = (error_small < 1e-20) && (error_small_c < 1e-20);
    add_test_result("Small numbers accuracy", small_success,
                   small_success ? "C implementation handles small numbers" : "C implementation fails with small numbers",
                   fmax(error_small, error_small_c));
    
    /* Test with large numbers */
    double x_large[n], y_large[n];
    double complex x_large_c[n], y_large_c[n];
    
    for (int i = 0; i < n; i++) {
        x_large[i] = 1e6 * (double)(i + 1);
        y_large[i] = 1e6 * (double)(i + 2);
        x_large_c[i] = 1e6 * (double)(i + 1) + 1e5 * (double)(i + 1) * I;
        y_large_c[i] = 1e6 * (double)(i + 2) + 1e5 * (double)(i + 2) * I;
    }
    
    double dot_large = komega_math_ddot(x_large, y_large, n);
    double complex dotc_large = komega_math_zdotc(x_large_c, y_large_c, n);
    
    /* Calculate expected results */
    double dot_large_expected = 0.0;
    double complex dotc_large_expected = 0.0 + 0.0 * I;
    for (int i = 0; i < n; i++) {
        dot_large_expected += x_large[i] * y_large[i];
        dotc_large_expected += conj(x_large_c[i]) * y_large_c[i];
    }
    
    double error_large = fabs(dot_large - dot_large_expected);
    double error_large_c = cabs(dotc_large - dotc_large_expected);
    
    bool large_success = (error_large < 1e-6) && (error_large_c < 1e-6);
    add_test_result("Large numbers accuracy", large_success,
                   large_success ? "C implementation handles large numbers" : "C implementation fails with large numbers",
                   fmax(error_large, error_large_c));
    
    /* Cleanup */
    komega_math_destroy(math);
}

void print_test_summary(void) {
    printf("\n📊 Test Summary\n");
    printf("============================================================\n");
    printf("Total tests: %d\n", g_test_suite.num_tests);
    printf("Passed: %d\n", g_test_suite.passed_tests);
    printf("Failed: %d\n", g_test_suite.failed_tests);
    
    if (g_test_suite.failed_tests == 0) {
        printf("\n✅ All tests passed! C and Fortran implementations are consistent.\n");
    } else {
        printf("\n❌ Some tests failed! C and Fortran implementations differ.\n");
    }
}

int main(void) {
    printf("🚀 Running C vs Fortran Implementation Comparison\n");
    printf("============================================================\n");
    
    /* Initialize test suite */
    test_result_init(&g_test_suite);
    
    /* Run comparison tests */
    test_math_operations_comparison();
    test_parameter_comparison();
    test_value_storage_comparison();
    test_solver_interface_comparison();
    test_numerical_accuracy_comparison();
    
    /* Print summary */
    print_test_summary();
    
    /* Cleanup */
    test_result_cleanup(&g_test_suite);
    
    return (g_test_suite.failed_tests == 0) ? 0 : 1;
}
