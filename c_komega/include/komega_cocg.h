/**
 * @file komega_cocg.h
 * @brief COCG solver header for Komega C library
 * 
 * This module provides the Conjugate Orthogonal CG (COCG) solver for complex
 * linear systems. It corresponds to the komega_cocg Fortran module.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 * 
 * @see komega_cocg.c
 * @see komega.h
 */

#ifndef KOMEGA_COCG_H
#define KOMEGA_COCG_H

#include "komega.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup COCGStruct COCG Solver Structure
 * @{
 */

/**
 * @brief COCG solver structure
 * 
 * This structure holds the state and parameters for the COCG solver.
 * It corresponds to the komega_cocg Fortran module.
 */
struct komega_cocg {
    bool initialized;           /**< Initialization flag */
    int ndim;                   /**< Problem dimension */
    int nl;                     /**< Local dimension */
    int nz;                     /**< Number of frequencies */
    int itermax;                /**< Maximum iterations */
    double threshold;           /**< Convergence threshold */
    int iter;                   /**< Current iteration */
    double resnorm;             /**< Residual norm */
    
    /* COCG-specific vectors */
    komega_complex_t* r;        /**< Residual vector */
    komega_complex_t* p;        /**< Search direction vector */
    komega_complex_t* v;        /**< Matrix-vector product vector */
    
    /* COCG coefficients */
    komega_complex_t alpha;     /**< Alpha coefficient */
    komega_complex_t beta;      /**< Beta coefficient */
    komega_complex_t rho;       /**< Rho coefficient */
    komega_complex_t rho_old;   /**< Previous rho coefficient */
    
    /* Frequency data */
    komega_complex_t* z;        /**< Frequency array */
    komega_complex_t z_seed;    /**< Seed frequency */
    
    /* Restart data */
    komega_complex_t* alpha_save; /**< Saved alpha coefficients */
    komega_complex_t* beta_save;  /**< Saved beta coefficients */
    komega_complex_t* r_l_save;   /**< Saved residual vectors */
};

/** @} */

/**
 * @defgroup COCGFunctions COCG Solver Functions
 * @{
 */

/**
 * @brief Create a new COCG solver instance
 * 
 * Allocates memory for a new COCG solver structure.
 * 
 * @return Pointer to new COCG solver instance, or NULL on failure
 * 
 * @see komega_cocg_destroy()
 */
struct komega_cocg* komega_cocg_create(void);

/**
 * @brief Destroy a COCG solver instance
 * 
 * Frees the memory allocated for the COCG solver structure.
 * 
 * @param solver Pointer to COCG solver instance to destroy
 * 
 * @see komega_cocg_create()
 */
void komega_cocg_destroy(struct komega_cocg* solver);

/**
 * @brief Initialize COCG solver
 * 
 * Initializes the COCG solver with the given parameters and vectors.
 * 
 * @param solver Pointer to COCG solver instance
 * @param ndim Problem dimension
 * @param nl Local dimension
 * @param nz Number of frequencies
 * @param x Solution vector (initialized to zero)
 * @param z Frequency array
 * @param itermax Maximum number of iterations
 * @param threshold Convergence threshold
 * @param comm MPI communicator (unused in single-threaded version)
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cocg_init(struct komega_cocg* solver,
                                 int ndim, int nl, int nz,
                                 komega_complex_t* x, komega_complex_t* z,
                                 int itermax, double threshold, int comm);

/**
 * @brief Restart COCG solver
 * 
 * Restarts the COCG solver with saved state information.
 * 
 * @param solver Pointer to COCG solver instance
 * @param ndim Problem dimension
 * @param nl Local dimension
 * @param nz Number of frequencies
 * @param x Solution vector
 * @param z Frequency array
 * @param itermax Maximum number of iterations
 * @param threshold Convergence threshold
 * @param status Solver status
 * @param iter_old Previous iteration count
 * @param v2 Vector v2
 * @param v12 Vector v12
 * @param alpha_save Saved alpha coefficients
 * @param beta_save Saved beta coefficients
 * @param z_seed Seed frequency
 * @param r_l_save Saved residual vectors
 * @param comm MPI communicator
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cocg_restart(struct komega_cocg* solver,
                                    int ndim, int nl, int nz,
                                    komega_complex_t* x, komega_complex_t* z,
                                    int itermax, double threshold,
                                    int* status, int* iter_old,
                                    komega_complex_t* v2, komega_complex_t* v12,
                                    komega_complex_t* alpha_save,
                                    komega_complex_t* beta_save,
                                    komega_complex_t* z_seed,
                                    komega_complex_t* r_l_save, int comm);

/**
 * @brief Update COCG solver iteration
 * 
 * Performs one iteration of the COCG solver.
 * 
 * @param solver Pointer to COCG solver instance
 * @param v12 Vector v12
 * @param v2 Vector v2
 * @param x Solution vector (updated)
 * @param r_l Residual vector (updated)
 * @param status Solver status (updated)
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cocg_update(struct komega_cocg* solver,
                                   komega_complex_t* v12, komega_complex_t* v2,
                                   komega_complex_t* x, komega_complex_t* r_l,
                                   int* status);

/**
 * @brief Get COCG solver coefficients
 * 
 * Retrieves the current solver coefficients.
 * 
 * @param solver Pointer to COCG solver instance
 * @param alpha_save Output array for alpha coefficients
 * @param beta_save Output array for beta coefficients
 * @param z_seed Output array for seed frequency
 * @param r_l_save Output array for saved residual vectors
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cocg_getcoef(struct komega_cocg* solver,
                                   komega_complex_t* alpha_save,
                                   komega_complex_t* beta_save,
                                   komega_complex_t* z_seed,
                                   komega_complex_t* r_l_save);

/**
 * @brief Get COCG solver vectors
 * 
 * Retrieves the current solver vectors.
 * 
 * @param solver Pointer to COCG solver instance
 * @param r_old Output array for old residual vector
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cocg_getvec(struct komega_cocg* solver, komega_complex_t* r_old);

/**
 * @brief Get COCG solver residual
 * 
 * Retrieves the current residual norm.
 * 
 * @param solver Pointer to COCG solver instance
 * @param res Output array for residual norms
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cocg_getresidual(struct komega_cocg* solver, double* res);

/**
 * @brief Finalize COCG solver
 * 
 * Cleans up the COCG solver and frees internal resources.
 * 
 * @param solver Pointer to COCG solver instance
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cocg_finalize(struct komega_cocg* solver);

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* KOMEGA_COCG_H */
