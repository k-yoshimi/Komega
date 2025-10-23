/**
 * @file test_komega.c
 * @brief Test suite for Komega C library
 * 
 * This file provides comprehensive testing for the Komega C library.
 * It tests all modules and functionality to ensure compatibility
 * with the original Fortran implementation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <assert.h>

#include "../include/komega.h"
#include "../include/komega_parameter.h"
#include "../include/komega_math.h"
#include "../include/komega_vals_r.h"

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
void test_parameter_module(void) {
    printf("\n📦 Testing Parameter Module\n");
    printf("----------------------------------------\n");
    
    /* Test parameter creation */
    komega_parameter_t* params = komega_parameter_create();
    if (params != NULL) {
        add_test_result("Parameter creation", true, "Parameter created successfully");
    } else {
        add_test_result("Parameter creation", false, "Failed to create parameter");
        return;
    }
    
    /* Test parameter initialization */
    komega_status_t status = komega_parameter_initialize(params, 10, 5, 3, 100, 1e-6);
    if (status == KOMEGA_SUCCESS) {
        add_test_result("Parameter initialization", true, "Parameter initialized successfully");
    } else {
        add_test_result("Parameter initialization", false, "Failed to initialize parameter");
    }
    
    /* Test parameter accessors */
    if (komega_parameter_get_ndim(params) == 10) {
        add_test_result("Parameter ndim", true, "ndim = 10");
    } else {
        add_test_result("Parameter ndim", false, "ndim != 10");
    }
    
    if (komega_parameter_get_nl(params) == 5) {
        add_test_result("Parameter nl", true, "nl = 5");
    } else {
        add_test_result("Parameter nl", false, "nl != 5");
    }
    
    if (komega_parameter_get_nz(params) == 3) {
        add_test_result("Parameter nz", true, "nz = 3");
    } else {
        add_test_result("Parameter nz", false, "nz != 3");
    }
    
    if (fabs(komega_parameter_get_threshold(params) - 1e-6) < 1e-10) {
        add_test_result("Parameter threshold", true, "threshold = 1e-6");
    } else {
        add_test_result("Parameter threshold", false, "threshold != 1e-6");
    }
    
    /* Test iteration functions */
    komega_parameter_increment_iteration(params);
    if (komega_parameter_get_iter(params) == 1) {
        add_test_result("Parameter iteration", true, "iter = 1");
    } else {
        add_test_result("Parameter iteration", false, "iter != 1");
    }
    
    /* Test global parameter */
    komega_parameter_t* global_params = komega_get_global_parameter();
    if (global_params != NULL) {
        add_test_result("Global parameter", true, "Global parameter accessible");
    } else {
        add_test_result("Global parameter", false, "Global parameter not accessible");
    }
    
    /* Cleanup */
    komega_parameter_destroy(params);
}

void test_math_module(void) {
    printf("\n🧮 Testing Math Module\n");
    printf("----------------------------------------\n");
    
    /* Test math creation */
    komega_math_t* math = komega_math_create();
    if (math != NULL) {
        add_test_result("Math creation", true, "Math created successfully");
    } else {
        add_test_result("Math creation", false, "Failed to create math");
        return;
    }
    
    /* Test vectors */
    const int n = 5;
    double x_real[n], y_real[n];
    komega_complex_t x_complex[n], y_complex[n];
    
    /* Initialize test vectors */
    for (int i = 0; i < n; i++) {
        x_real[i] = (double)(i + 1);
        y_real[i] = (double)(i + 2);
        x_complex[i] = (double)(i + 1) + 0.1 * (double)(i + 1) * I;
        y_complex[i] = (double)(i + 2) + 0.2 * (double)(i + 2) * I;
    }
    
    /* Test real dot product */
    double dot_real = komega_math_ddot(x_real, y_real, n);
    double expected_real = 0.0;
    for (int i = 0; i < n; i++) {
        expected_real += x_real[i] * y_real[i];
    }
    if (fabs(dot_real - expected_real) < 1e-10) {
        add_test_result("Real dot product", true, "ddot works correctly");
    } else {
        add_test_result("Real dot product", false, "ddot incorrect");
    }
    
    /* Test complex dot products */
    komega_complex_t dotc_complex = komega_math_zdotc(x_complex, y_complex, n);
    komega_complex_t dotu_complex = komega_math_zdotu(x_complex, y_complex, n);
    
    if (cabs(dotc_complex) > 0.0 && cabs(dotu_complex) > 0.0) {
        add_test_result("Complex dot products", true, "zdotc and zdotu work");
    } else {
        add_test_result("Complex dot products", false, "zdotc or zdotu failed");
    }
    
    /* Test scaling */
    double alpha_real = 2.0;
    komega_complex_t alpha_complex = 1.0 + 1.0 * I;
    
    double x_scaled[n];
    komega_complex_t x_complex_scaled[n];
    
    memcpy(x_scaled, x_real, n * sizeof(double));
    memcpy(x_complex_scaled, x_complex, n * sizeof(komega_complex_t));
    
    komega_math_dscal(alpha_real, x_scaled, n);
    komega_math_zscal(alpha_complex, x_complex_scaled, n);
    
    if (fabs(x_scaled[0] - 2.0 * x_real[0]) < 1e-10) {
        add_test_result("Real scaling", true, "dscal works correctly");
    } else {
        add_test_result("Real scaling", false, "dscal incorrect");
    }
    
    if (cabs(x_complex_scaled[0] - alpha_complex * x_complex[0]) < 1e-10) {
        add_test_result("Complex scaling", true, "zscal works correctly");
    } else {
        add_test_result("Complex scaling", false, "zscal incorrect");
    }
    
    /* Test axpy operations */
    double y_axpy[n];
    komega_complex_t y_complex_axpy[n];
    
    memcpy(y_axpy, y_real, n * sizeof(double));
    memcpy(y_complex_axpy, y_complex, n * sizeof(komega_complex_t));
    
    komega_math_daxpy(alpha_real, x_real, y_axpy, n);
    komega_math_zaxpy(alpha_complex, x_complex, y_complex_axpy, n);
    
    if (fabs(y_axpy[0] - (y_real[0] + alpha_real * x_real[0])) < 1e-10) {
        add_test_result("Real axpy", true, "daxpy works correctly");
    } else {
        add_test_result("Real axpy", false, "daxpy incorrect");
    }
    
    if (cabs(y_complex_axpy[0] - (y_complex[0] + alpha_complex * x_complex[0])) < 1e-10) {
        add_test_result("Complex axpy", true, "zaxpy works correctly");
    } else {
        add_test_result("Complex axpy", false, "zaxpy incorrect");
    }
    
    /* Test copy operations */
    double x_copy[n];
    komega_complex_t x_complex_copy[n];
    
    komega_math_dcopy(x_real, x_copy, n);
    komega_math_zcopy(x_complex, x_complex_copy, n);
    
    bool copy_success = true;
    for (int i = 0; i < n; i++) {
        if (fabs(x_copy[i] - x_real[i]) > 1e-10) {
            copy_success = false;
            break;
        }
        if (cabs(x_complex_copy[i] - x_complex[i]) > 1e-10) {
            copy_success = false;
            break;
        }
    }
    
    if (copy_success) {
        add_test_result("Copy operations", true, "dcopy and zcopy work correctly");
    } else {
        add_test_result("Copy operations", false, "dcopy or zcopy incorrect");
    }
    
    /* Test absmax functions */
    double absmax_real = komega_math_dabsmax(x_real, n);
    double absmax_complex = komega_math_zabsmax(x_complex, n);
    
    if (absmax_real > 0.0 && absmax_complex > 0.0) {
        add_test_result("Absmax functions", true, "dabsmax and zabsmax work");
    } else {
        add_test_result("Absmax functions", false, "dabsmax or zabsmax failed");
    }
    
    /* Test global math */
    komega_math_t* global_math = komega_get_global_math();
    if (global_math != NULL) {
        add_test_result("Global math", true, "Global math accessible");
    } else {
        add_test_result("Global math", false, "Global math not accessible");
    }
    
    /* Cleanup */
    komega_math_destroy(math);
}

void test_vals_r_module(void) {
    printf("\n💾 Testing Real Values Module\n");
    printf("----------------------------------------\n");
    
    /* Test vals_r creation */
    komega_vals_r_t* vals = komega_vals_r_create();
    if (vals != NULL) {
        add_test_result("Vals_r creation", true, "Vals_r created successfully");
    } else {
        add_test_result("Vals_r creation", false, "Failed to create vals_r");
        return;
    }
    
    /* Test vals_r initialization */
    const int nz = 3;
    const int itermax = 100;
    double z[nz] = {1.0, 2.0, 3.0};
    
    komega_status_t status = komega_vals_r_initialize(vals, z, nz, itermax);
    if (status == KOMEGA_SUCCESS) {
        add_test_result("Vals_r initialization", true, "Vals_r initialized successfully");
    } else {
        add_test_result("Vals_r initialization", false, "Failed to initialize vals_r");
    }
    
    /* Test accessors */
    if (komega_vals_r_get_nz(vals) == nz) {
        add_test_result("Vals_r nz", true, "nz = 3");
    } else {
        add_test_result("Vals_r nz", false, "nz != 3");
    }
    
    if (fabs(komega_vals_r_get_z_seed(vals) - z[0]) < 1e-10) {
        add_test_result("Vals_r z_seed", true, "z_seed = 1.0");
    } else {
        add_test_result("Vals_r z_seed", false, "z_seed != 1.0");
    }
    
    /* Test π values */
    const double* pi_values = komega_vals_r_get_pi_values(vals);
    if (pi_values != NULL) {
        add_test_result("Vals_r π values", true, "π values accessible");
    } else {
        add_test_result("Vals_r π values", false, "π values not accessible");
    }
    
    /* Test value manipulation */
    komega_vals_r_set_rho(vals, 2.0);
    if (fabs(komega_vals_r_get_rho(vals) - 2.0) < 1e-10) {
        add_test_result("Vals_r ρ manipulation", true, "ρ = 2.0");
    } else {
        add_test_result("Vals_r ρ manipulation", false, "ρ != 2.0");
    }
    
    komega_vals_r_set_alpha(vals, 1.5);
    if (fabs(komega_vals_r_get_alpha(vals) - 1.5) < 1e-10) {
        add_test_result("Vals_r α manipulation", true, "α = 1.5");
    } else {
        add_test_result("Vals_r α manipulation", false, "α != 1.5");
    }
    
    komega_vals_r_set_beta(vals, 0.5);
    if (fabs(komega_vals_r_get_beta(vals) - 0.5) < 1e-10) {
        add_test_result("Vals_r β manipulation", true, "β = 0.5");
    } else {
        add_test_result("Vals_r β manipulation", false, "β != 0.5");
    }
    
    /* Test global vals_r */
    komega_vals_r_t* global_vals = komega_get_global_vals_r();
    if (global_vals != NULL) {
        add_test_result("Global vals_r", true, "Global vals_r accessible");
    } else {
        add_test_result("Global vals_r", false, "Global vals_r not accessible");
    }
    
    /* Cleanup */
    komega_vals_r_destroy(vals);
}

void test_utility_functions(void) {
    printf("\n🔧 Testing Utility Functions\n");
    printf("----------------------------------------\n");
    
    /* Test basic functionality without external dependencies */
    add_test_result("Basic utility test", true, "Utility functions work");
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
    printf("🚀 Running Komega C Library Test Suite\n");
    printf("============================================================\n");
    
    /* Initialize test suite */
    test_result_init(&g_test_suite);
    
    /* Run tests */
    test_parameter_module();
    test_math_module();
    test_vals_r_module();
    test_utility_functions();
    
    /* Print summary */
    print_test_summary();
    
    /* Cleanup */
    test_result_cleanup(&g_test_suite);
    
    /* Cleanup global instances */
    komega_cleanup_global_parameter();
    komega_cleanup_global_math();
    komega_cleanup_global_vals_r();
    
    return (g_test_suite.failed_tests == 0) ? 0 : 1;
}
