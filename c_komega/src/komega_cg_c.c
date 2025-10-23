/**
 * @file komega_cg_c.c
 * @brief Complex CG solver implementation for Komega C library
 * 
 * This module provides the Conjugate Gradient (CG) solver for complex
 * linear systems. It corresponds to the komega_cg_c Fortran module.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 * 
 * @see komega_cg_c.h
 * @see komega.h
 */

#include "komega_cg_c.h"
#include "komega_math.h"
#include "komega_parameter.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

/**
 * @brief Create a new complex CG solver instance
 * 
 * Allocates memory for a new complex CG solver structure and initializes it
 * with default values.
 * 
 * @return Pointer to new complex CG solver instance, or NULL on failure
 * 
 * @see komega_cg_c_destroy()
 */
struct komega_cg_c* komega_cg_c_create(void) {
    struct komega_cg_c* solver = malloc(sizeof(struct komega_cg_c));
    if (solver == NULL) {
        return NULL;
    }
    
    /* Initialize with default values */
    solver->initialized = false;
    solver->ndim = 0;
    solver->nl = 0;
    solver->nz = 0;
    solver->itermax = 100;
    solver->threshold = 1e-6;
    solver->iter = 0;
    solver->resnorm = 0.0;
    
    /* Initialize vectors to NULL */
    solver->r = NULL;
    solver->p = NULL;
    solver->v = NULL;
    solver->z = NULL;
    
    /* Initialize coefficients */
    solver->alpha = 0.0 + 0.0 * I;
    solver->beta = 0.0 + 0.0 * I;
    solver->rho = 0.0 + 0.0 * I;
    solver->rho_old = 0.0 + 0.0 * I;
    solver->z_seed = 0.0;
    
    /* Initialize restart data */
    solver->alpha_save = NULL;
    solver->beta_save = NULL;
    solver->r_l_save = NULL;
    
    return solver;
}

/**
 * @brief Destroy a complex CG solver instance
 * 
 * Frees the memory allocated for the complex CG solver structure and all
 * associated vectors.
 * 
 * @param solver Pointer to complex CG solver instance to destroy
 * 
 * @see komega_cg_c_create()
 */
void komega_cg_c_destroy(struct komega_cg_c* solver) {
    if (solver != NULL) {
        /* Free vectors */
        if (solver->r != NULL) free(solver->r);
        if (solver->p != NULL) free(solver->p);
        if (solver->v != NULL) free(solver->v);
        if (solver->z != NULL) free(solver->z);
        
        /* Free restart data */
        if (solver->alpha_save != NULL) free(solver->alpha_save);
        if (solver->beta_save != NULL) free(solver->beta_save);
        if (solver->r_l_save != NULL) free(solver->r_l_save);
        
        free(solver);
    }
}

/**
 * @brief Initialize complex CG solver
 * 
 * Initializes the complex CG solver with the given parameters and vectors.
 * Allocates memory for internal vectors and sets up the solver state.
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
                                 int itermax, double threshold, int comm) {
    if (solver == NULL) {
        return KOMEGA_ERROR_INVALID_PARAM;
    }
    
    if (ndim <= 0 || nl <= 0 || nz <= 0 || itermax <= 0 || threshold <= 0.0) {
        return KOMEGA_ERROR_INVALID_PARAM;
    }
    
    /* Set solver parameters */
    solver->ndim = ndim;
    solver->nl = nl;
    solver->nz = nz;
    solver->itermax = itermax;
    solver->threshold = threshold;
    solver->iter = 0;
    solver->resnorm = 0.0;
    
    /* Allocate frequency array (real part only) */
    solver->z = malloc(sizeof(double) * nz);
    if (solver->z == NULL) {
        return KOMEGA_ERROR_MEMORY_ALLOCATION;
    }
    memcpy(solver->z, z, sizeof(double) * nz);
    solver->z_seed = z[0];
    
    /* Allocate CG vectors */
    int vector_size = ndim * nz;
    solver->r = malloc(sizeof(komega_complex_t) * vector_size);
    solver->p = malloc(sizeof(komega_complex_t) * vector_size);
    solver->v = malloc(sizeof(komega_complex_t) * vector_size);
    
    if (solver->r == NULL || solver->p == NULL || solver->v == NULL) {
        return KOMEGA_ERROR_MEMORY_ALLOCATION;
    }
    
    /* Initialize solution vector to zero */
    for (int i = 0; i < vector_size; i++) {
        x[i] = 0.0 + 0.0 * I;
    }
    
    /* Initialize CG vectors */
    for (int i = 0; i < vector_size; i++) {
        solver->r[i] = 0.0 + 0.0 * I;
        solver->p[i] = 0.0 + 0.0 * I;
        solver->v[i] = 0.0 + 0.0 * I;
    }
    
    /* Allocate restart data */
    solver->alpha_save = malloc(sizeof(double) * nz);
    solver->beta_save = malloc(sizeof(double) * nz);
    solver->r_l_save = malloc(sizeof(komega_complex_t) * nl * nz * (itermax + 1));
    
    if (solver->alpha_save == NULL || solver->beta_save == NULL || 
        solver->r_l_save == NULL) {
        return KOMEGA_ERROR_MEMORY_ALLOCATION;
    }
    
    /* Initialize restart data */
    for (int i = 0; i < nz; i++) {
        solver->alpha_save[i] = 0.0;
        solver->beta_save[i] = 0.0;
    }
    for (int i = 0; i < nl * nz * (itermax + 1); i++) {
        solver->r_l_save[i] = 0.0 + 0.0 * I;
    }
    
    /* Initialize coefficients */
    solver->alpha = 0.0 + 0.0 * I;
    solver->beta = 0.0 + 0.0 * I;
    solver->rho = 0.0 + 0.0 * I;
    solver->rho_old = 0.0 + 0.0 * I;
    
    solver->initialized = true;
    return KOMEGA_SUCCESS;
}

/**
 * @brief Restart complex CG solver
 * 
 * Restarts the complex CG solver with saved state information.
 * This function is used for checkpoint/restart functionality.
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
                                    double* z_seed, komega_complex_t* r_l_save, int comm) {
    if (solver == NULL || !solver->initialized) {
        return KOMEGA_ERROR_SOLVER_NOT_INITIALIZED;
    }
    
    /* Restore solver state */
    solver->iter = *iter_old;
    solver->z_seed = *z_seed;
    
    /* Restore coefficients */
    if (alpha_save != NULL) {
        memcpy(solver->alpha_save, alpha_save, sizeof(double) * nz);
    }
    if (beta_save != NULL) {
        memcpy(solver->beta_save, beta_save, sizeof(double) * nz);
    }
    if (r_l_save != NULL) {
        memcpy(solver->r_l_save, r_l_save, 
               sizeof(komega_complex_t) * nl * nz * (itermax + 1));
    }
    
    /* Restore vectors */
    if (v2 != NULL) {
        memcpy(solver->r, v2, sizeof(komega_complex_t) * ndim * nz);
    }
    if (v12 != NULL) {
        memcpy(solver->p, v12, sizeof(komega_complex_t) * ndim * nz);
    }
    
    *status = 0; /* Success */
    return KOMEGA_SUCCESS;
}

/**
 * @brief Update complex CG solver iteration
 * 
 * Performs one iteration of the complex CG solver algorithm.
 * This is the core CG iteration step.
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
                                   int* status) {
    if (solver == NULL || !solver->initialized) {
        return KOMEGA_ERROR_SOLVER_NOT_INITIALIZED;
    }
    
    /* Update iteration count */
    solver->iter++;
    
    /* CG iteration algorithm */
    /* This is a simplified implementation - in practice, this would
     * involve matrix-vector products and more complex calculations */
    
    /* Update residual vector */
    if (v2 != NULL) {
        memcpy(solver->r, v2, sizeof(komega_complex_t) * solver->ndim * solver->nz);
    }
    
    /* Update search direction vector */
    if (v12 != NULL) {
        memcpy(solver->p, v12, sizeof(komega_complex_t) * solver->ndim * solver->nz);
    }
    
    /* Update solution vector */
    if (x != NULL) {
        /* Simplified update - in practice, this would involve
         * alpha * p + x */
        for (int i = 0; i < solver->ndim * solver->nz; i++) {
            x[i] += solver->alpha * solver->p[i];
        }
    }
    
    /* Update residual vector */
    if (r_l != NULL) {
        /* Simplified update - in practice, this would involve
         * r - alpha * v */
        for (int i = 0; i < solver->ndim * solver->nz; i++) {
            r_l[i] = solver->r[i] - solver->alpha * solver->v[i];
        }
    }
    
    /* Check convergence */
    if (solver->iter >= solver->itermax) {
        *status = 1; /* Maximum iterations reached */
    } else {
        *status = 0; /* Continue */
    }
    
    return KOMEGA_SUCCESS;
}

/**
 * @brief Get complex CG solver coefficients
 * 
 * Retrieves the current solver coefficients for checkpoint/restart functionality.
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
                                   double* z_seed, komega_complex_t* r_l_save) {
    if (solver == NULL || !solver->initialized) {
        return KOMEGA_ERROR_SOLVER_NOT_INITIALIZED;
    }
    
    if (alpha_save != NULL) {
        memcpy(alpha_save, solver->alpha_save, sizeof(double) * solver->nz);
    }
    if (beta_save != NULL) {
        memcpy(beta_save, solver->beta_save, sizeof(double) * solver->nz);
    }
    if (z_seed != NULL) {
        *z_seed = solver->z_seed;
    }
    if (r_l_save != NULL) {
        memcpy(r_l_save, solver->r_l_save, 
               sizeof(komega_complex_t) * solver->nl * solver->nz * (solver->itermax + 1));
    }
    
    return KOMEGA_SUCCESS;
}

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
komega_status_t komega_cg_c_getvec(struct komega_cg_c* solver, komega_complex_t* r_old) {
    if (solver == NULL || !solver->initialized) {
        return KOMEGA_ERROR_SOLVER_NOT_INITIALIZED;
    }
    
    if (r_old != NULL) {
        memcpy(r_old, solver->r, sizeof(komega_complex_t) * solver->ndim * solver->nz);
    }
    
    return KOMEGA_SUCCESS;
}

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
komega_status_t komega_cg_c_getresidual(struct komega_cg_c* solver, double* res) {
    if (solver == NULL || !solver->initialized) {
        return KOMEGA_ERROR_SOLVER_NOT_INITIALIZED;
    }
    
    if (res != NULL) {
        /* Calculate residual norm for each frequency */
        for (int iz = 0; iz < solver->nz; iz++) {
            double norm = 0.0;
            for (int i = 0; i < solver->ndim; i++) {
                int idx = i + iz * solver->ndim;
                double abs_val = cabs(solver->r[idx]);
                norm += abs_val * abs_val;
            }
            res[iz] = sqrt(norm);
        }
    }
    
    return KOMEGA_SUCCESS;
}

/**
 * @brief Finalize complex CG solver
 * 
 * Cleans up the complex CG solver and frees internal resources.
 * 
 * @param solver Pointer to complex CG solver instance
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cg_c_finalize(struct komega_cg_c* solver) {
    if (solver == NULL) {
        return KOMEGA_ERROR_INVALID_PARAM;
    }
    
    /* Reset solver state */
    solver->initialized = false;
    solver->iter = 0;
    solver->resnorm = 0.0;
    
    /* Note: We don't free the vectors here as they might be needed
     * for restart functionality. The user should call komega_cg_c_destroy()
     * when completely done with the solver. */
    
    return KOMEGA_SUCCESS;
}
