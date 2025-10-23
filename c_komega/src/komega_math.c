/**
 * @file komega_math.c
 * @brief Mathematical operations module for Komega C library
 * 
 * This module provides mathematical operations similar to BLAS functions.
 * It corresponds to the komega_math Fortran module and includes
 * dot products, scaling, axpy operations, and other mathematical utilities.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 * 
 * @see komega_math.h
 * @see komega.h
 */

#include "komega_math.h"
#include <math.h>
#include <stdlib.h>

/* Global math instance */
static komega_math_t* global_math = NULL;

/**
 * @brief Create a new math instance
 * 
 * Allocates memory for a new math structure and initializes it.
 * 
 * @return Pointer to new math instance, or NULL on failure
 * 
 * @see komega_math_destroy()
 */
komega_math_t* komega_math_create(void) {
    komega_math_t* math = malloc(sizeof(komega_math_t));
    if (math == NULL) {
        return NULL;
    }
    
    /* Initialize with default values */
    math->initialized = true;
    
    return math;
}

/**
 * @brief Destroy a math instance
 * 
 * Frees the memory allocated for the math structure.
 * Safe to call with NULL pointer.
 * 
 * @param math Pointer to math instance to destroy
 * 
 * @see komega_math_create()
 */
void komega_math_destroy(komega_math_t* math) {
    if (math != NULL) {
        free(math);
    }
}

/**
 * @brief Real dot product
 * 
 * Computes the dot product of two real vectors: result = x^T * y
 * 
 * @param x First vector
 * @param y Second vector
 * @param n Vector length
 * @return Dot product result
 */
double komega_math_ddot(const double* x, const double* y, int n) {
    double result = 0.0;
    for (int i = 0; i < n; i++) {
        result += x[i] * y[i];
    }
    return result;
}

/**
 * @brief Complex conjugate dot product
 * 
 * Computes the conjugate dot product of two complex vectors: result = x^H * y
 * 
 * @param x First vector
 * @param y Second vector
 * @param n Vector length
 * @return Conjugate dot product result
 */
komega_complex_t komega_math_zdotc(const komega_complex_t* x, 
                                   const komega_complex_t* y, int n) {
    komega_complex_t result = 0.0 + 0.0 * I;
    for (int i = 0; i < n; i++) {
        result += conj(x[i]) * y[i];
    }
    return result;
}

/**
 * @brief Complex dot product
 * 
 * Computes the dot product of two complex vectors: result = x^T * y
 * 
 * @param x First vector
 * @param y Second vector
 * @param n Vector length
 * @return Dot product result
 */
komega_complex_t komega_math_zdotu(const komega_complex_t* x, 
                                   const komega_complex_t* y, int n) {
    komega_complex_t result = 0.0 + 0.0 * I;
    for (int i = 0; i < n; i++) {
        result += x[i] * y[i];
    }
    return result;
}

/**
 * @brief Real vector scaling
 * 
 * Scales a real vector by a scalar: x = alpha * x
 * 
 * @param alpha Scaling factor
 * @param x Vector to scale (modified in place)
 * @param n Vector length
 */
void komega_math_dscal(double alpha, double* x, int n) {
    for (int i = 0; i < n; i++) {
        x[i] *= alpha;
    }
}

/**
 * @brief Complex vector scaling
 * 
 * Scales a complex vector by a scalar: x = alpha * x
 * 
 * @param alpha Scaling factor
 * @param x Vector to scale (modified in place)
 * @param n Vector length
 */
void komega_math_zscal(komega_complex_t alpha, komega_complex_t* x, int n) {
    for (int i = 0; i < n; i++) {
        x[i] *= alpha;
    }
}

/**
 * @brief Real axpy operation
 * 
 * Performs y = alpha * x + y
 * 
 * @param alpha Scaling factor
 * @param x First vector
 * @param y Second vector (modified in place)
 * @param n Vector length
 */
void komega_math_daxpy(double alpha, const double* x, double* y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] += alpha * x[i];
    }
}

/**
 * @brief Complex axpy operation
 * 
 * Performs y = alpha * x + y
 * 
 * @param alpha Scaling factor
 * @param x First vector
 * @param y Second vector (modified in place)
 * @param n Vector length
 */
void komega_math_zaxpy(komega_complex_t alpha, const komega_complex_t* x, 
                       komega_complex_t* y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] += alpha * x[i];
    }
}

/**
 * @brief Real vector copy
 * 
 * Copies a real vector: y = x
 * 
 * @param x Source vector
 * @param y Destination vector
 * @param n Vector length
 */
void komega_math_dcopy(const double* x, double* y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = x[i];
    }
}

/**
 * @brief Complex vector copy
 * 
 * Copies a complex vector: y = x
 * 
 * @param x Source vector
 * @param y Destination vector
 * @param n Vector length
 */
void komega_math_zcopy(const komega_complex_t* x, komega_complex_t* y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = x[i];
    }
}

/**
 * @brief Real vector absolute maximum
 * 
 * Finds the maximum absolute value in a real vector
 * 
 * @param x Vector
 * @param n Vector length
 * @return Maximum absolute value
 */
double komega_math_dabsmax(const double* x, int n) {
    double max_val = 0.0;
    for (int i = 0; i < n; i++) {
        double abs_val = fabs(x[i]);
        if (abs_val > max_val) {
            max_val = abs_val;
        }
    }
    return max_val;
}

/**
 * @brief Complex vector absolute maximum
 * 
 * Finds the maximum absolute value in a complex vector
 * 
 * @param x Vector
 * @param n Vector length
 * @return Maximum absolute value
 */
double komega_math_zabsmax(const komega_complex_t* x, int n) {
    double max_val = 0.0;
    for (int i = 0; i < n; i++) {
        double abs_val = cabs(x[i]);
        if (abs_val > max_val) {
            max_val = abs_val;
        }
    }
    return max_val;
}

/**
 * @brief Real dot product (MPI-like interface)
 * 
 * MPI-like interface for real dot product computation.
 * Currently simulated for single-threaded execution.
 * 
 * @param x First vector
 * @param y Second vector
 * @param n Vector length
 * @return Dot product result
 */
double komega_math_ddot_mpi(const double* x, const double* y, int n) {
    return komega_math_ddot(x, y, n);
}

/**
 * @brief Complex conjugate dot product (MPI-like interface)
 * 
 * MPI-like interface for complex conjugate dot product computation.
 * Currently simulated for single-threaded execution.
 * 
 * @param x First vector
 * @param y Second vector
 * @param n Vector length
 * @return Conjugate dot product result
 */
komega_complex_t komega_math_zdotc_mpi(const komega_complex_t* x, 
                                       const komega_complex_t* y, int n) {
    return komega_math_zdotc(x, y, n);
}

/**
 * @brief Complex dot product (MPI-like interface)
 * 
 * MPI-like interface for complex dot product computation.
 * Currently simulated for single-threaded execution.
 * 
 * @param x First vector
 * @param y Second vector
 * @param n Vector length
 * @return Dot product result
 */
komega_complex_t komega_math_zdotu_mpi(const komega_complex_t* x, 
                                       const komega_complex_t* y, int n) {
    return komega_math_zdotu(x, y, n);
}

/**
 * @brief Get the global math instance
 * 
 * Returns a singleton instance of the math structure.
 * Creates a new instance if none exists.
 * 
 * @return Pointer to global math instance
 * 
 * @see komega_cleanup_global_math()
 */
komega_math_t* komega_get_global_math(void) {
    if (global_math == NULL) {
        global_math = komega_math_create();
    }
    return global_math;
}

/**
 * @brief Clean up the global math instance
 * 
 * Destroys the global math instance and frees its memory.
 * 
 * @see komega_get_global_math()
 */
void komega_cleanup_global_math(void) {
    if (global_math != NULL) {
        komega_math_destroy(global_math);
        global_math = NULL;
    }
}
