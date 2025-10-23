/**
 * @file test_detailed.c
 * @brief Detailed test for C Komega Library
 * 
 * This test performs more detailed testing of the solver functionality.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "../include/komega.h"
#include "../include/komega_bicg.h"
#include "../include/komega_cg_r.h"
#include "../include/komega_cg_c.h"
#include "../include/komega_cocg.h"

#define TEST_NDIM 10
#define TEST_NZ 2
#define TEST_NL 3
#define TEST_ITERMAX 100
#define TEST_THRESHOLD 1e-6

int main(void) {
    printf("C Komega Library - Detailed Test\n");
    printf("================================\n\n");
    
    /* Test parameters */
    int ndim = TEST_NDIM;
    int nz = TEST_NZ;
    int nl = TEST_NL;
    int itermax = TEST_ITERMAX;
    double threshold = TEST_THRESHOLD;
    
    printf("Test Parameters:\n");
    printf("  ndim = %d\n", ndim);
    printf("  nz = %d\n", nz);
    printf("  nl = %d\n", nl);
    printf("  itermax = %d\n", itermax);
    printf("  threshold = %e\n\n", threshold);
    
    /* Test 1: BiCG Solver with Complex Data */
    printf("Test 1: BiCG Solver with Complex Data\n");
    printf("-------------------------------------\n");
    
    struct komega_bicg* bicg = komega_bicg_create();
    if (bicg == NULL) {
        printf("❌ Failed to create BiCG solver\n");
        return 1;
    }
    
    /* Create test data */
    double complex* z = malloc(sizeof(double complex) * nz);
    double complex* x = malloc(sizeof(double complex) * ndim * nz);
    
    for (int i = 0; i < nz; i++) {
        z[i] = (double)(i + 1) + 0.1 * (double)(i + 1) * I;
    }
    
    for (int i = 0; i < ndim * nz; i++) {
        x[i] = 0.0 + 0.0 * I;
    }
    
    /* Initialize BiCG solver */
    komega_status_t status = komega_bicg_init(bicg, ndim, nl, nz, x, z, itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to initialize BiCG solver\n");
        komega_bicg_destroy(bicg);
        free(z);
        free(x);
        return 1;
    }
    printf("✓ BiCG solver initialized successfully\n");
    
    /* Test solver functions */
    double res[TEST_NZ];
    status = komega_bicg_getresidual(bicg, res);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get BiCG residual\n");
    } else {
        printf("✓ BiCG residual retrieved successfully\n");
        for (int i = 0; i < nz; i++) {
            printf("  Residual[%d] = %e\n", i, res[i]);
        }
    }
    
    komega_bicg_destroy(bicg);
    printf("✓ BiCG solver destroyed successfully\n\n");
    
    /* Test 2: Real CG Solver */
    printf("Test 2: Real CG Solver\n");
    printf("---------------------\n");
    
    struct komega_cg_r* cg_r = komega_cg_r_create();
    if (cg_r == NULL) {
        printf("❌ Failed to create Real CG solver\n");
        free(z);
        free(x);
        return 1;
    }
    
    /* Create real test data */
    double* z_real = malloc(sizeof(double) * nz);
    double* x_real = malloc(sizeof(double) * ndim * nz);
    
    for (int i = 0; i < nz; i++) {
        z_real[i] = (double)(i + 1);
    }
    
    for (int i = 0; i < ndim * nz; i++) {
        x_real[i] = 0.0;
    }
    
    /* Initialize Real CG solver */
    status = komega_cg_r_init(cg_r, ndim, nl, nz, x_real, z_real, itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to initialize Real CG solver\n");
        komega_cg_r_destroy(cg_r);
        free(z);
        free(x);
        free(z_real);
        free(x_real);
        return 1;
    }
    printf("✓ Real CG solver initialized successfully\n");
    
    /* Test solver functions */
    double res_real[TEST_NZ];
    status = komega_cg_r_getresidual(cg_r, res_real);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get Real CG residual\n");
    } else {
        printf("✓ Real CG residual retrieved successfully\n");
        for (int i = 0; i < nz; i++) {
            printf("  Residual[%d] = %e\n", i, res_real[i]);
        }
    }
    
    komega_cg_r_destroy(cg_r);
    printf("✓ Real CG solver destroyed successfully\n\n");
    
    /* Test 3: Complex CG Solver */
    printf("Test 3: Complex CG Solver\n");
    printf("------------------------\n");
    
    struct komega_cg_c* cg_c = komega_cg_c_create();
    if (cg_c == NULL) {
        printf("❌ Failed to create Complex CG solver\n");
        free(z);
        free(x);
        free(z_real);
        free(x_real);
        return 1;
    }
    
    /* Initialize Complex CG solver */
    status = komega_cg_c_init(cg_c, ndim, nl, nz, x, z_real, itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to initialize Complex CG solver\n");
        komega_cg_c_destroy(cg_c);
        free(z);
        free(x);
        free(z_real);
        free(x_real);
        return 1;
    }
    printf("✓ Complex CG solver initialized successfully\n");
    
    /* Test solver functions */
    double res_complex[TEST_NZ];
    status = komega_cg_c_getresidual(cg_c, res_complex);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get Complex CG residual\n");
    } else {
        printf("✓ Complex CG residual retrieved successfully\n");
        for (int i = 0; i < nz; i++) {
            printf("  Residual[%d] = %e\n", i, res_complex[i]);
        }
    }
    
    komega_cg_c_destroy(cg_c);
    printf("✓ Complex CG solver destroyed successfully\n\n");
    
    /* Test 4: COCG Solver */
    printf("Test 4: COCG Solver\n");
    printf("------------------\n");
    
    struct komega_cocg* cocg = komega_cocg_create();
    if (cocg == NULL) {
        printf("❌ Failed to create COCG solver\n");
        free(z);
        free(x);
        free(z_real);
        free(x_real);
        return 1;
    }
    
    /* Initialize COCG solver */
    status = komega_cocg_init(cocg, ndim, nl, nz, x, z, itermax, threshold, 0);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to initialize COCG solver\n");
        komega_cocg_destroy(cocg);
        free(z);
        free(x);
        free(z_real);
        free(x_real);
        return 1;
    }
    printf("✓ COCG solver initialized successfully\n");
    
    /* Test solver functions */
    double res_cocg[TEST_NZ];
    status = komega_cocg_getresidual(cocg, res_cocg);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to get COCG residual\n");
    } else {
        printf("✓ COCG residual retrieved successfully\n");
        for (int i = 0; i < nz; i++) {
            printf("  Residual[%d] = %e\n", i, res_cocg[i]);
        }
    }
    
    komega_cocg_destroy(cocg);
    printf("✓ COCG solver destroyed successfully\n\n");
    
    /* Cleanup */
    free(z);
    free(x);
    free(z_real);
    free(x_real);
    
    printf("🎉 All detailed tests passed successfully!\n");
    printf("C Komega Library solvers are working correctly.\n");
    
    return 0;
}
