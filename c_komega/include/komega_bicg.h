/**
 * @file komega_bicg.h
 * @brief BiCG solver header for Komega C library
 * 
 * This module provides the Bi-Conjugate Gradient (BiCG) solver for complex
 * linear systems. It corresponds to the komega_bicg Fortran module.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 * 
 * @see komega_bicg.c
 * @see komega.h
 */

#ifndef KOMEGA_BICG_H
#define KOMEGA_BICG_H

#include "komega.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup BiCGStruct BiCG Solver Structure
 * @{
 */

/**
 * @brief BiCG solver structure
 * 
 * This structure holds the state and parameters for the BiCG solver.
 * It corresponds to the komega_bicg Fortran module.
 */
struct komega_bicg {
    bool initialized;           /**< Initialization flag */
    int ndim;                   /**< Problem dimension */
    int nl;                     /**< Local dimension */
    int nz;                     /**< Number of frequencies */
    int itermax;                /**< Maximum iterations */
    double threshold;           /**< Convergence threshold */
    int iter;                   /**< Current iteration */
    double resnorm;             /**< Residual norm */
    
    /* BiCG-specific vectors */
    komega_complex_t* r;        /**< Residual vector */
    komega_complex_t* r_tilde;  /**< Bi-orthogonal residual vector */
    komega_complex_t* p;        /**< Search direction vector */
    komega_complex_t* p_tilde;  /**< Bi-orthogonal search direction vector */
    komega_complex_t* v;        /**< Matrix-vector product vector */
    komega_complex_t* v_tilde;  /**< Bi-orthogonal matrix-vector product vector */
    
    /* BiCG coefficients */
    komega_complex_t alpha;     /**< Alpha coefficient */
    komega_complex_t beta;      /**< Beta coefficient */
    komega_complex_t rho;       /**< Rho coefficient */
    komega_complex_t rho_old;   /**< Previous rho coefficient */
    komega_complex_t omega;     /**< Omega coefficient */
    komega_complex_t omega_old; /**< Previous omega coefficient */
    
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
 * @defgroup BiCGFunctions BiCG Solver Functions
 * @{
 */

/**
 * @brief Create a new BiCG solver instance
 * 
 * Allocates memory for a new BiCG solver structure.
 * 
 * @return Pointer to new BiCG solver instance, or NULL on failure
 * 
 * @see komega_bicg_destroy()
 */
struct komega_bicg* komega_bicg_create(void);

/**
 * @brief Destroy a BiCG solver instance
 * 
 * Frees the memory allocated for the BiCG solver structure.
 * 
 * @param solver Pointer to BiCG solver instance to destroy
 * 
 * @see komega_bicg_create()
 */
void komega_bicg_destroy(struct komega_bicg* solver);

/**
 * @brief Initialize BiCG solver
 * 
 * Initializes the BiCG solver with the given parameters and vectors.
 * 
 * @param solver Pointer to BiCG solver instance
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
komega_status_t komega_bicg_init(struct komega_bicg* solver,
                                int ndim, int nl, int nz,
                                komega_complex_t* x, komega_complex_t* z,
                                int itermax, double threshold, int comm);

/**
 * @brief Restart BiCG solver
 * 
 * Restarts the BiCG solver with saved state information.
 * 
 * @param solver Pointer to BiCG solver instance
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
 * @param v4 Vector v4
 * @param v14 Vector v14
 * @param alpha_save Saved alpha coefficients
 * @param beta_save Saved beta coefficients
 * @param z_seed Seed frequency
 * @param r_l_save Saved residual vectors
 * @param comm MPI communicator
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_bicg_restart(struct komega_bicg* solver,
                                    int ndim, int nl, int nz,
                                    komega_complex_t* x, komega_complex_t* z,
                                    int itermax, double threshold,
                                    int* status, int* iter_old,
                                    komega_complex_t* v2, komega_complex_t* v12,
                                    komega_complex_t* v4, komega_complex_t* v14,
                                    komega_complex_t* alpha_save,
                                    komega_complex_t* beta_save,
                                    komega_complex_t* z_seed,
                                    komega_complex_t* r_l_save, int comm);

/**
 * @brief Update BiCG solver iteration
 * 
 * Performs one iteration of the BiCG solver.
 * 
 * @param solver Pointer to BiCG solver instance
 * @param v12 Vector v12
 * @param v2 Vector v2
 * @param v14 Vector v14
 * @param v4 Vector v4
 * @param x Solution vector (updated)
 * @param r_l Residual vector (updated)
 * @param status Solver status (updated)
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_bicg_update(struct komega_bicg* solver,
                                   komega_complex_t* v12, komega_complex_t* v2,
                                   komega_complex_t* v14, komega_complex_t* v4,
                                   komega_complex_t* x, komega_complex_t* r_l,
                                   int* status);

/**
 * @brief Get BiCG solver coefficients
 * 
 * Retrieves the current solver coefficients.
 * 
 * @param solver Pointer to BiCG solver instance
 * @param alpha_save Output array for alpha coefficients
 * @param beta_save Output array for beta coefficients
 * @param z_seed Output array for seed frequency
 * @param r_l_save Output array for saved residual vectors
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_bicg_getcoef(struct komega_bicg* solver,
                                   komega_complex_t* alpha_save,
                                   komega_complex_t* beta_save,
                                   komega_complex_t* z_seed,
                                   komega_complex_t* r_l_save);

/**
 * @brief Get BiCG solver vectors
 * 
 * Retrieves the current solver vectors.
 * 
 * @param solver Pointer to BiCG solver instance
 * @param r_old Output array for old residual vector
 * @param r_tilde_old Output array for old bi-orthogonal residual vector
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_bicg_getvec(struct komega_bicg* solver,
                                  komega_complex_t* r_old,
                                  komega_complex_t* r_tilde_old);

/**
 * @brief Get BiCG solver residual
 * 
 * Retrieves the current residual norm.
 * 
 * @param solver Pointer to BiCG solver instance
 * @param res Output array for residual norms
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_bicg_getresidual(struct komega_bicg* solver, double* res);

/**
 * @brief Finalize BiCG solver
 * 
 * Cleans up the BiCG solver and frees internal resources.
 * 
 * @param solver Pointer to BiCG solver instance
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_bicg_finalize(struct komega_bicg* solver);

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* KOMEGA_BICG_H */
