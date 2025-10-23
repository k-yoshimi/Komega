/**
 * @file komega_parameter.c
 * @brief Parameter module for Komega C library
 * 
 * This module provides global parameter management for the Komega library.
 * It corresponds to the komega_parameter Fortran module and handles
 * solver parameters, iteration counts, and convergence thresholds.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 * 
 * @see komega_parameter.h
 * @see komega.h
 */

#include "komega_parameter.h"
#include <stdlib.h>
#include <string.h>

/* Global parameter instance */
static komega_parameter_t* global_params = NULL;

/**
 * @brief Create a new parameter instance
 * 
 * Allocates memory for a new parameter structure and initializes it
 * with default values. The default values are:
 * - ndim = 0
 * - nl = 0  
 * - nz = 0
 * - itermax = 100
 * - threshold = 1e-6
 * - almost0 = KOMEGA_ALMOST_ZERO
 * - iter = 0
 * - resnorm = 0.0
 * - initialized = false
 * 
 * @return Pointer to new parameter instance, or NULL on failure
 * 
 * @see komega_parameter_destroy()
 * @see komega_parameter_initialize()
 */
komega_parameter_t* komega_parameter_create(void) {
    komega_parameter_t* params = malloc(sizeof(komega_parameter_t));
    if (params == NULL) {
        return NULL;
    }
    
    /* Initialize with default values */
    params->ndim = 0;
    params->nl = 0;
    params->nz = 0;
    params->itermax = 100;
    params->threshold = 1e-6;
    params->almost0 = KOMEGA_ALMOST_ZERO;
    params->iter = 0;
    params->resnorm = 0.0;
    params->initialized = false;
    
    return params;
}

/**
 * @brief Destroy a parameter instance
 * 
 * Frees the memory allocated for the parameter structure.
 * Safe to call with NULL pointer.
 * 
 * @param params Pointer to parameter instance to destroy
 * 
 * @see komega_parameter_create()
 */
void komega_parameter_destroy(komega_parameter_t* params) {
    if (params != NULL) {
        free(params);
    }
}

/**
 * @brief Initialize parameter instance with solver parameters
 * 
 * Sets up the parameter structure with the given solver parameters.
 * All parameters are validated before setting. The function will
 * return an error if any parameter is invalid (<= 0).
 * 
 * @param params Pointer to parameter instance
 * @param ndim Dimension of the problem (must be > 0)
 * @param nl Local dimension (must be > 0)
 * @param nz Number of frequencies (must be > 0)
 * @param itermax Maximum number of iterations (must be > 0)
 * @param threshold Convergence threshold (must be > 0)
 * 
 * @return KOMEGA_SUCCESS on success, KOMEGA_ERROR_INVALID_PARAM on failure
 * 
 * @see komega_parameter_create()
 */
komega_status_t komega_parameter_initialize(komega_parameter_t* params, 
                                           int ndim, int nl, int nz, 
                                           int itermax, double threshold) {
    if (params == NULL) {
        return KOMEGA_ERROR_INVALID_PARAM;
    }
    
    if (ndim <= 0 || nl <= 0 || nz <= 0 || itermax <= 0 || threshold <= 0.0) {
        return KOMEGA_ERROR_INVALID_PARAM;
    }
    
    params->ndim = ndim;
    params->nl = nl;
    params->nz = nz;
    params->itermax = itermax;
    params->threshold = threshold;
    params->iter = 0;
    params->resnorm = 0.0;
    params->initialized = true;
    
    return KOMEGA_SUCCESS;
}

/**
 * @brief Increment the iteration counter
 * 
 * Increments the current iteration count by 1.
 * Safe to call with NULL pointer.
 * 
 * @param params Pointer to parameter instance
 */
void komega_parameter_increment_iteration(komega_parameter_t* params) {
    if (params != NULL) {
        params->iter++;
    }
}

/**
 * @brief Reset the iteration counter
 * 
 * Resets the current iteration count to 0.
 * Safe to call with NULL pointer.
 * 
 * @param params Pointer to parameter instance
 */
void komega_parameter_reset_iteration(komega_parameter_t* params) {
    if (params != NULL) {
        params->iter = 0;
    }
}

/**
 * @brief Check if parameters are initialized
 * 
 * @param params Pointer to parameter instance
 * @return true if initialized, false otherwise
 */
bool komega_parameter_is_initialized(const komega_parameter_t* params) {
    return (params != NULL) ? params->initialized : false;
}

/**
 * @brief Get the problem dimension
 * 
 * @param params Pointer to parameter instance
 * @return Problem dimension (ndim), or 0 if params is NULL
 */
int komega_parameter_get_ndim(const komega_parameter_t* params) {
    return (params != NULL) ? params->ndim : 0;
}

/**
 * @brief Get the local dimension
 * 
 * @param params Pointer to parameter instance
 * @return Local dimension (nl), or 0 if params is NULL
 */
int komega_parameter_get_nl(const komega_parameter_t* params) {
    return (params != NULL) ? params->nl : 0;
}

/**
 * @brief Get the number of frequencies
 * 
 * @param params Pointer to parameter instance
 * @return Number of frequencies (nz), or 0 if params is NULL
 */
int komega_parameter_get_nz(const komega_parameter_t* params) {
    return (params != NULL) ? params->nz : 0;
}

/**
 * @brief Get the maximum number of iterations
 * 
 * @param params Pointer to parameter instance
 * @return Maximum number of iterations (itermax), or 0 if params is NULL
 */
int komega_parameter_get_itermax(const komega_parameter_t* params) {
    return (params != NULL) ? params->itermax : 0;
}

/**
 * @brief Get the convergence threshold
 * 
 * @param params Pointer to parameter instance
 * @return Convergence threshold, or 0.0 if params is NULL
 */
double komega_parameter_get_threshold(const komega_parameter_t* params) {
    return (params != NULL) ? params->threshold : 0.0;
}

/**
 * @brief Get the current iteration count
 * 
 * @param params Pointer to parameter instance
 * @return Current iteration count, or 0 if params is NULL
 */
int komega_parameter_get_iter(const komega_parameter_t* params) {
    return (params != NULL) ? params->iter : 0;
}

/**
 * @brief Get the residual norm
 * 
 * @param params Pointer to parameter instance
 * @return Current residual norm, or 0.0 if params is NULL
 */
double komega_parameter_get_resnorm(const komega_parameter_t* params) {
    return (params != NULL) ? params->resnorm : 0.0;
}

/**
 * @brief Set the residual norm
 * 
 * @param params Pointer to parameter instance
 * @param resnorm New residual norm value
 */
void komega_parameter_set_resnorm(komega_parameter_t* params, double resnorm) {
    if (params != NULL) {
        params->resnorm = resnorm;
    }
}

/**
 * @brief Get the global parameter instance
 * 
 * Returns a singleton instance of the parameter structure.
 * Creates a new instance if none exists.
 * 
 * @return Pointer to global parameter instance
 * 
 * @see komega_cleanup_global_parameter()
 */
komega_parameter_t* komega_get_global_parameter(void) {
    if (global_params == NULL) {
        global_params = komega_parameter_create();
    }
    return global_params;
}

/**
 * @brief Clean up the global parameter instance
 * 
 * Destroys the global parameter instance and frees its memory.
 * 
 * @see komega_get_global_parameter()
 */
void komega_cleanup_global_parameter(void) {
    if (global_params != NULL) {
        komega_parameter_destroy(global_params);
        global_params = NULL;
    }
}
