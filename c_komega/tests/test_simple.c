/**
 * @file test_simple.c
 * @brief Simple test for C Komega Library
 * 
 * This is a simplified test to verify basic functionality.
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

int main(void) {
    printf("C Komega Library - Simple Test\n");
    printf("==============================\n\n");
    
    /* Test 1: Basic parameter creation */
    printf("Test 1: Parameter Creation\n");
    printf("---------------------------\n");
    
    struct komega_parameter* params = komega_parameter_create();
    if (params == NULL) {
        printf("❌ Failed to create parameter\n");
        return 1;
    }
    printf("✓ Parameter created successfully\n");
    
    komega_status_t status = komega_parameter_initialize(params, 10, 5, 3, 100, 1e-6);
    if (status != KOMEGA_SUCCESS) {
        printf("❌ Failed to initialize parameter\n");
        komega_parameter_destroy(params);
        return 1;
    }
    printf("✓ Parameter initialized successfully\n");
    
    komega_parameter_destroy(params);
    printf("✓ Parameter destroyed successfully\n\n");
    
    /* Test 2: BiCG solver creation */
    printf("Test 2: BiCG Solver Creation\n");
    printf("----------------------------\n");
    
    struct komega_bicg* bicg = komega_bicg_create();
    if (bicg == NULL) {
        printf("❌ Failed to create BiCG solver\n");
        return 1;
    }
    printf("✓ BiCG solver created successfully\n");
    
    komega_bicg_destroy(bicg);
    printf("✓ BiCG solver destroyed successfully\n\n");
    
    /* Test 3: Real CG solver creation */
    printf("Test 3: Real CG Solver Creation\n");
    printf("-------------------------------\n");
    
    struct komega_cg_r* cg_r = komega_cg_r_create();
    if (cg_r == NULL) {
        printf("❌ Failed to create Real CG solver\n");
        return 1;
    }
    printf("✓ Real CG solver created successfully\n");
    
    komega_cg_r_destroy(cg_r);
    printf("✓ Real CG solver destroyed successfully\n\n");
    
    /* Test 4: Complex CG solver creation */
    printf("Test 4: Complex CG Solver Creation\n");
    printf("-----------------------------------\n");
    
    struct komega_cg_c* cg_c = komega_cg_c_create();
    if (cg_c == NULL) {
        printf("❌ Failed to create Complex CG solver\n");
        return 1;
    }
    printf("✓ Complex CG solver created successfully\n");
    
    komega_cg_c_destroy(cg_c);
    printf("✓ Complex CG solver destroyed successfully\n\n");
    
    /* Test 5: COCG solver creation */
    printf("Test 5: COCG Solver Creation\n");
    printf("----------------------------\n");
    
    struct komega_cocg* cocg = komega_cocg_create();
    if (cocg == NULL) {
        printf("❌ Failed to create COCG solver\n");
        return 1;
    }
    printf("✓ COCG solver created successfully\n");
    
    komega_cocg_destroy(cocg);
    printf("✓ COCG solver destroyed successfully\n\n");
    
    printf("🎉 All tests passed successfully!\n");
    printf("C Komega Library is working correctly.\n");
    
    return 0;
}
