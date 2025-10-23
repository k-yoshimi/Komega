/**
 * @file komega_cg_c.h
 * @brief Complex CG solver header for Komega C library
 * 
 * This module provides the Conjugate Gradient (CG) solver for complex
 * linear systems. It corresponds to the komega_cg_c Fortran module.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 * 
 * @see komega_cg_c.c
 * @see komega.h
 */

#ifndef KOMEGA_CG_C_H
#define KOMEGA_CG_C_H

#include "komega.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup CGCStruct Complex CG Solver Structure
 * @{
 */

/**
 * @brief Complex CG solver structure
 * 
 * This structure holds the state and parameters for the complex CG solver.
 * It corresponds to the komega_cg_c Fortran module.
 */
struct komega_cg_c {
    bool initialized;           /**< Initialization flag */
    int ndim;                   /**< Problem dimension */
    int nl;                     /**< Local dimension */
    int nz;                     /**< Number of frequencies */
    int itermax;                /**< Maximum iterations */
    double threshold;           /**< Convergence threshold */
    int iter;                   /**< Current iteration */
    double resnorm;             /**< Residual norm */
    
    /* CG-specific vectors */
    komega_complex_t* r;        /**< Residual vector */
    komega_complex_t* p;        /**< Search direction vector */
    komega_complex_t* v;        /**< Matrix-vector product vector */
    
    /* CG coefficients */
    komega_complex_t alpha;     /**< Alpha coefficient */
    komega_complex_t beta;      /**< Beta coefficient */
    komega_complex_t rho;       /**< Rho coefficient */
    komega_complex_t rho_old;   /**< Previous rho coefficient */
    
    /* Frequency data */
    double* z;                  /**< Frequency array (real part only) */
    double z_seed;              /**< Seed frequency (real part) */
    
    /* Restart data */
    double* alpha_save;         /**< Saved alpha coefficients (real) */
    double* beta_save;          /**< Saved beta coefficients (real) */
    komega_complex_t* r_l_save; /**< Saved residual vectors */
};

/** @} */

/**
 * @defgroup CGCFunctions Complex CG Solver Functions
 * @{
 */

/**
 * @brief Create a new complex CG solver instance
 * 
 * Allocates memory for a new complex CG solver structure.
 * 
 * @return Pointer to new complex CG solver instance, or NULL on failure
 * 
 * @see komega_cg_c_destroy()
 */
struct komega_cg_c* komega_cg_c_create(void);

/**
 * @brief Destroy a complex CG solver instance
 * 
 * Frees the memory allocated for the complex CG solver structure.
 * 
 * @param solver Pointer to complex CG solver instance to destroy
 * 
 * @see komega_cg_c_create()
 */
void komega_cg_c_destroy(struct komega_cg_c* solver);

/**
 * @brief Initialize complex CG solver
 * 
 * Initializes the complex CG solver with the given parameters and vectors.
 * 
 * @param solver Pointer to complex CG solver instance
 * @param ndim Problem dimension
 * @param nl Local dimension
 * @param nz Number of frequencies
 * @param x Solution vector (initialized to zero)
 * @param z Frequency array (real part only)
 * @param itermax Maximum number of iterations
 * @param threshold Convergence threshold
 * @param comm MPI communicator (unused in single-threaded version)
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_c_init(struct komega_cg_c* solver,
                                 int ndim, int nl, int nz,
                                 komega_complex_t* x, double* z,
                                 int itermax, double threshold, int comm);

/**
 * @brief Restart complex CG solver
 * 
 * Restarts the complex CG solver with saved state information.
 * 
 * @param solver Pointer to complex CG solver instance
 * @param ndim Problem dimension
 * @param nl Local dimension
 * @param nz Number of frequencies
 * @param x Solution vector
 * @param z Frequency array (real part only)
 * @param itermax Maximum number of iterations
 * @param threshold Convergence threshold
 * @param status Solver status
 * @param iter_old Previous iteration count
 * @param v2 Vector v2
 * @param v12 Vector v12
 * @param alpha_save Saved alpha coefficients (real)
 * @param beta_save Saved beta coefficients (real)
 * @param z_seed Seed frequency (real)
 * @param r_l_save Saved residual vectors
 * @param comm MPI communicator
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_c_restart(struct komega_cg_c* solver,
                                    int ndim, int nl, int nz,
                                    komega_complex_t* x, double* z,
                                    int itermax, double threshold,
                                    int* status, int* iter_old,
                                    komega_complex_t* v2, komega_complex_t* v12,
                                    double* alpha_save, double* beta_save,
                                    double* z_seed, komega_complex_t* r_l_save, int comm);

/**
 * @brief Update complex CG solver iteration
 * 
 * Performs one iteration of the complex CG solver.
 * 
 * @param solver Pointer to complex CG solver instance
 * @param v12 Vector v12
 * @param v2 Vector v2
 * @param x Solution vector (updated)
 * @param r_l Residual vector (updated)
 * @param status Solver status (updated)
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_c_update(struct komega_cg_c* solver,
                                   komega_complex_t* v12, komega_complex_t* v2,
                                   komega_complex_t* x, komega_complex_t* r_l,
                                   int* status);

/**
 * @brief Get complex CG solver coefficients
 * 
 * Retrieves the current solver coefficients.
 * 
 * @param solver Pointer to complex CG solver instance
 * @param alpha_save Output array for alpha coefficients (real)
 * @param beta_save Output array for beta coefficients (real)
 * @param z_seed Output array for seed frequency (real)
 * @param r_l_save Output array for saved residual vectors
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_c_getcoef(struct komega_cg_c* solver,
                                   double* alpha_save, double* beta_save,
                                   double* z_seed, komega_complex_t* r_l_save);

/**
 * @brief Get complex CG solver vectors
 * 
 * Retrieves the current solver vectors.
 * 
 * @param solver Pointer to complex CG solver instance
 * @param r_old Output array for old residual vector
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_c_getvec(struct komega_cg_c* solver, komega_complex_t* r_old);

/**
 * @brief Get complex CG solver residual
 * 
 * Retrieves the current residual norm.
 * 
 * @param solver Pointer to complex CG solver instance
 * @param res Output array for residual norms
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_c_getresidual(struct komega_cg_c* solver, double* res);

/**
 * @brief Finalize complex CG solver
 * 
 * Cleans up the complex CG solver and frees internal resources.
 * 
 * @param solver Pointer to complex CG solver instance
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_c_finalize(struct komega_cg_c* solver);

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* KOMEGA_CG_C_H */
