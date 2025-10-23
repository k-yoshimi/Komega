/**
 * @file komega_parameter.h
 * @brief Parameter module header for Komega C library
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
 * @see komega_parameter.c
 * @see komega.h
 */

#ifndef KOMEGA_PARAMETER_H
#define KOMEGA_PARAMETER_H

#include "komega.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup ParameterStruct Parameter Structure
 * @{
 */

/**
 * @brief Parameter structure for Komega library
 * 
 * This structure holds all the global parameters needed for the Komega solvers.
 * It corresponds to the komega_parameter Fortran module.
 */
struct komega_parameter {
    int ndim;           /**< Dimension of the problem */
    int nl;             /**< Local dimension */
    int nz;             /**< Number of frequencies */
    int itermax;        /**< Maximum number of iterations */
    double threshold;   /**< Convergence threshold */
    double almost0;     /**< Almost zero threshold */
    int iter;           /**< Current iteration */
    double resnorm;     /**< Residual norm */
    bool initialized;   /**< Initialization flag */
};

/** @} */

/**
 * @defgroup ParameterFunctions Parameter Functions
 * @{
 */

/**
 * @brief Create a new parameter instance
 * 
 * Allocates memory for a new parameter structure and initializes it
 * with default values.
 * 
 * @return Pointer to new parameter instance, or NULL on failure
 * 
 * @see komega_parameter_destroy()
 */
komega_parameter_t* komega_parameter_create(void);

/**
 * @brief Destroy a parameter instance
 * 
 * Frees the memory allocated for the parameter structure.
 * 
 * @param params Pointer to parameter instance to destroy
 * 
 * @see komega_parameter_create()
 */
void komega_parameter_destroy(komega_parameter_t* params);

/**
 * @brief Initialize parameter instance with solver parameters
 * 
 * Sets up the parameter structure with the given solver parameters.
 * All parameters are validated before setting.
 * 
 * @param params Pointer to parameter instance
 * @param ndim Dimension of the problem (must be > 0)
 * @param nl Local dimension (must be > 0)
 * @param nz Number of frequencies (must be > 0)
 * @param itermax Maximum number of iterations (must be > 0)
 * @param threshold Convergence threshold (must be > 0)
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 * 
 * @see komega_parameter_create()
 */
komega_status_t komega_parameter_initialize(komega_parameter_t* params, 
                                           int ndim, int nl, int nz, 
                                           int itermax, double threshold);

/**
 * @defgroup ParameterAccessors Parameter Accessor Functions
 * @{
 */

/**
 * @brief Increment the iteration counter
 * 
 * Increments the current iteration count by 1.
 * 
 * @param params Pointer to parameter instance
 */
void komega_parameter_increment_iteration(komega_parameter_t* params);

/**
 * @brief Reset the iteration counter
 * 
 * Resets the current iteration count to 0.
 * 
 * @param params Pointer to parameter instance
 */
void komega_parameter_reset_iteration(komega_parameter_t* params);

/**
 * @brief Check if parameters are initialized
 * 
 * @param params Pointer to parameter instance
 * @return true if initialized, false otherwise
 */
bool komega_parameter_is_initialized(const komega_parameter_t* params);

/**
 * @brief Get the problem dimension
 * 
 * @param params Pointer to parameter instance
 * @return Problem dimension (ndim)
 */
int komega_parameter_get_ndim(const komega_parameter_t* params);

/**
 * @brief Get the local dimension
 * 
 * @param params Pointer to parameter instance
 * @return Local dimension (nl)
 */
int komega_parameter_get_nl(const komega_parameter_t* params);

/**
 * @brief Get the number of frequencies
 * 
 * @param params Pointer to parameter instance
 * @return Number of frequencies (nz)
 */
int komega_parameter_get_nz(const komega_parameter_t* params);

/**
 * @brief Get the maximum number of iterations
 * 
 * @param params Pointer to parameter instance
 * @return Maximum number of iterations (itermax)
 */
int komega_parameter_get_itermax(const komega_parameter_t* params);

/**
 * @brief Get the convergence threshold
 * 
 * @param params Pointer to parameter instance
 * @return Convergence threshold
 */
double komega_parameter_get_threshold(const komega_parameter_t* params);

/**
 * @brief Get the current iteration count
 * 
 * @param params Pointer to parameter instance
 * @return Current iteration count
 */
int komega_parameter_get_iter(const komega_parameter_t* params);

/**
 * @brief Get the residual norm
 * 
 * @param params Pointer to parameter instance
 * @return Current residual norm
 */
double komega_parameter_get_resnorm(const komega_parameter_t* params);

/**
 * @brief Set the residual norm
 * 
 * @param params Pointer to parameter instance
 * @param resnorm New residual norm value
 */
void komega_parameter_set_resnorm(komega_parameter_t* params, double resnorm);

/** @} */

/**
 * @defgroup ParameterGlobal Global Parameter Functions
 * @{
 */

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
komega_parameter_t* komega_get_global_parameter(void);

/**
 * @brief Clean up the global parameter instance
 * 
 * Destroys the global parameter instance and frees its memory.
 * 
 * @see komega_get_global_parameter()
 */
void komega_cleanup_global_parameter(void);

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* KOMEGA_PARAMETER_H */
