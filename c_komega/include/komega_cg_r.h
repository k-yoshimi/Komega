/**
 * @file komega_cg_r.h
 * @brief Real CG solver header for Komega C library
 * 
 * This module provides the Conjugate Gradient (CG) solver for real
 * linear systems. It corresponds to the komega_cg_r Fortran module.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 * 
 * @see komega_cg_r.c
 * @see komega.h
 */

#ifndef KOMEGA_CG_R_H
#define KOMEGA_CG_R_H

#include "komega.h"

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup CGRStruct Real CG Solver Structure
 * @{
 */

/**
 * @brief Real CG solver structure
 * 
 * This structure holds the state and parameters for the real CG solver.
 * It corresponds to the komega_cg_r Fortran module.
 */
struct komega_cg_r {
    bool initialized;           /**< Initialization flag */
    int ndim;                   /**< Problem dimension */
    int nl;                     /**< Local dimension */
    int nz;                     /**< Number of frequencies */
    int itermax;                /**< Maximum iterations */
    double threshold;           /**< Convergence threshold */
    int iter;                   /**< Current iteration */
    double resnorm;             /**< Residual norm */
    
    /* CG-specific vectors */
    double* r;                  /**< Residual vector */
    double* p;                  /**< Search direction vector */
    double* v;                  /**< Matrix-vector product vector */
    
    /* CG coefficients */
    double alpha;               /**< Alpha coefficient */
    double beta;                /**< Beta coefficient */
    double rho;                 /**< Rho coefficient */
    double rho_old;             /**< Previous rho coefficient */
    
    /* Frequency data */
    double* z;                  /**< Frequency array (real) */
    double z_seed;              /**< Seed frequency */
    
    /* Restart data */
    double* alpha_save;         /**< Saved alpha coefficients */
    double* beta_save;         /**< Saved beta coefficients */
    double* r_l_save;          /**< Saved residual vectors */
};

/** @} */

/**
 * @defgroup CGRFunctions Real CG Solver Functions
 * @{
 */

/**
 * @brief Create a new real CG solver instance
 * 
 * Allocates memory for a new real CG solver structure.
 * 
 * @return Pointer to new real CG solver instance, or NULL on failure
 * 
 * @see komega_cg_r_destroy()
 */
struct komega_cg_r* komega_cg_r_create(void);

/**
 * @brief Destroy a real CG solver instance
 * 
 * Frees the memory allocated for the real CG solver structure.
 * 
 * @param solver Pointer to real CG solver instance to destroy
 * 
 * @see komega_cg_r_create()
 */
void komega_cg_r_destroy(struct komega_cg_r* solver);

/**
 * @brief Initialize real CG solver
 * 
 * Initializes the real CG solver with the given parameters and vectors.
 * 
 * @param solver Pointer to real CG solver instance
 * @param ndim Problem dimension
 * @param nl Local dimension
 * @param nz Number of frequencies
 * @param x Solution vector (initialized to zero)
 * @param z Frequency array (real)
 * @param itermax Maximum number of iterations
 * @param threshold Convergence threshold
 * @param comm MPI communicator (unused in single-threaded version)
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_r_init(struct komega_cg_r* solver,
                                 int ndim, int nl, int nz,
                                 double* x, double* z,
                                 int itermax, double threshold, int comm);

/**
 * @brief Restart real CG solver
 * 
 * Restarts the real CG solver with saved state information.
 * 
 * @param solver Pointer to real CG solver instance
 * @param ndim Problem dimension
 * @param nl Local dimension
 * @param nz Number of frequencies
 * @param x Solution vector
 * @param z Frequency array (real)
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
komega_status_t komega_cg_r_restart(struct komega_cg_r* solver,
                                    int ndim, int nl, int nz,
                                    double* x, double* z,
                                    int itermax, double threshold,
                                    int* status, int* iter_old,
                                    double* v2, double* v12,
                                    double* alpha_save, double* beta_save,
                                    double* z_seed, double* r_l_save, int comm);

/**
 * @brief Update real CG solver iteration
 * 
 * Performs one iteration of the real CG solver.
 * 
 * @param solver Pointer to real CG solver instance
 * @param v12 Vector v12
 * @param v2 Vector v2
 * @param x Solution vector (updated)
 * @param r_l Residual vector (updated)
 * @param status Solver status (updated)
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_r_update(struct komega_cg_r* solver,
                                   double* v12, double* v2,
                                   double* x, double* r_l, int* status);

/**
 * @brief Get real CG solver coefficients
 * 
 * Retrieves the current solver coefficients.
 * 
 * @param solver Pointer to real CG solver instance
 * @param alpha_save Output array for alpha coefficients
 * @param beta_save Output array for beta coefficients
 * @param z_seed Output array for seed frequency
 * @param r_l_save Output array for saved residual vectors
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_r_getcoef(struct komega_cg_r* solver,
                                   double* alpha_save, double* beta_save,
                                   double* z_seed, double* r_l_save);

/**
 * @brief Get real CG solver vectors
 * 
 * Retrieves the current solver vectors.
 * 
 * @param solver Pointer to real CG solver instance
 * @param r_old Output array for old residual vector
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_r_getvec(struct komega_cg_r* solver, double* r_old);

/**
 * @brief Get real CG solver residual
 * 
 * Retrieves the current residual norm.
 * 
 * @param solver Pointer to real CG solver instance
 * @param res Output array for residual norms
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_r_getresidual(struct komega_cg_r* solver, double* res);

/**
 * @brief Finalize real CG solver
 * 
 * Cleans up the real CG solver and frees internal resources.
 * 
 * @param solver Pointer to real CG solver instance
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_r_finalize(struct komega_cg_r* solver);

/** @} */

#ifdef __cplusplus
}
#endif

#endif /* KOMEGA_CG_R_H */
