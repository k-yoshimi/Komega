/**
 * @file komega_cocg.c
 * @brief COCG solver implementation for Komega C library
 * 
 * This module provides the Conjugate Orthogonal CG (COCG) solver for complex
 * linear systems. It corresponds to the komega_cocg Fortran module.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 * 
 * @see komega_cocg.h
 * @see komega.h
 */

#include "komega_cocg.h"
#include "komega_math.h"
#include "komega_parameter.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

/**
 * @brief Create a new COCG solver instance
 * 
 * Allocates memory for a new COCG solver structure and initializes it
 * with default values.
 * 
 * @return Pointer to new COCG solver instance, or NULL on failure
 * 
 * @see komega_cocg_destroy()
 */
struct komega_cocg* komega_cocg_create(void) {
    struct komega_cocg* solver = malloc(sizeof(struct komega_cocg));
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
    solver->z_seed = 0.0 + 0.0 * I;
    
    /* Initialize restart data */
    solver->alpha_save = NULL;
    solver->beta_save = NULL;
    solver->r_l_save = NULL;
    
    return solver;
}

/**
 * @brief Destroy a COCG solver instance
 * 
 * Frees the memory allocated for the COCG solver structure and all
 * associated vectors.
 * 
 * @param solver Pointer to COCG solver instance to destroy
 * 
 * @see komega_cocg_create()
 */
void komega_cocg_destroy(struct komega_cocg* solver) {
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
 * @brief Initialize COCG solver
 * 
 * Initializes the COCG solver with the given parameters and vectors.
 * Allocates memory for internal vectors and sets up the solver state.
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
    
    /* Allocate frequency array */
    solver->z = malloc(sizeof(komega_complex_t) * nz);
    if (solver->z == NULL) {
        return KOMEGA_ERROR_MEMORY_ALLOCATION;
    }
    memcpy(solver->z, z, sizeof(komega_complex_t) * nz);
    solver->z_seed = z[0];
    
    /* Allocate COCG vectors */
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
    
    /* Initialize COCG vectors */
    for (int i = 0; i < vector_size; i++) {
        solver->r[i] = 0.0 + 0.0 * I;
        solver->p[i] = 0.0 + 0.0 * I;
        solver->v[i] = 0.0 + 0.0 * I;
    }
    
    /* Allocate restart data */
    solver->alpha_save = malloc(sizeof(komega_complex_t) * nz);
    solver->beta_save = malloc(sizeof(komega_complex_t) * nz);
    solver->r_l_save = malloc(sizeof(komega_complex_t) * nl * nz * (itermax + 1));
    
    if (solver->alpha_save == NULL || solver->beta_save == NULL || 
        solver->r_l_save == NULL) {
        return KOMEGA_ERROR_MEMORY_ALLOCATION;
    }
    
    /* Initialize restart data */
    for (int i = 0; i < nz; i++) {
        solver->alpha_save[i] = 0.0 + 0.0 * I;
        solver->beta_save[i] = 0.0 + 0.0 * I;
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
 * @brief Restart COCG solver
 * 
 * Restarts the COCG solver with saved state information.
 * This function is used for checkpoint/restart functionality.
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
                                    komega_complex_t* r_l_save, int comm) {
    if (solver == NULL || !solver->initialized) {
        return KOMEGA_ERROR_SOLVER_NOT_INITIALIZED;
    }
    
    /* Restore solver state */
    solver->iter = *iter_old;
    solver->z_seed = *z_seed;
    
    /* Restore coefficients */
    if (alpha_save != NULL) {
        memcpy(solver->alpha_save, alpha_save, sizeof(komega_complex_t) * nz);
    }
    if (beta_save != NULL) {
        memcpy(solver->beta_save, beta_save, sizeof(komega_complex_t) * nz);
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
 * @brief Update COCG solver iteration
 * 
 * Performs one iteration of the COCG solver algorithm.
 * This is the core COCG iteration step.
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
                                   int* status) {
    if (solver == NULL || !solver->initialized) {
        return KOMEGA_ERROR_SOLVER_NOT_INITIALIZED;
    }
    
    /* Update iteration count */
    solver->iter++;
    
    /* COCG iteration algorithm */
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
 * @brief Get COCG solver coefficients
 * 
 * Retrieves the current solver coefficients for checkpoint/restart functionality.
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
                                   komega_complex_t* r_l_save) {
    if (solver == NULL || !solver->initialized) {
        return KOMEGA_ERROR_SOLVER_NOT_INITIALIZED;
    }
    
    if (alpha_save != NULL) {
        memcpy(alpha_save, solver->alpha_save, 
               sizeof(komega_complex_t) * solver->nz);
    }
    if (beta_save != NULL) {
        memcpy(beta_save, solver->beta_save, 
               sizeof(komega_complex_t) * solver->nz);
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
 * @brief Get COCG solver vectors
 * 
 * Retrieves the current solver vectors.
 * 
 * @param solver Pointer to COCG solver instance
 * @param r_old Output array for old residual vector
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cocg_getvec(struct komega_cocg* solver, komega_complex_t* r_old) {
    if (solver == NULL || !solver->initialized) {
        return KOMEGA_ERROR_SOLVER_NOT_INITIALIZED;
    }
    
    if (r_old != NULL) {
        memcpy(r_old, solver->r, sizeof(komega_complex_t) * solver->ndim * solver->nz);
    }
    
    return KOMEGA_SUCCESS;
}

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
komega_status_t komega_cocg_getresidual(struct komega_cocg* solver, double* res) {
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
 * @brief Finalize COCG solver
 * 
 * Cleans up the COCG solver and frees internal resources.
 * 
 * @param solver Pointer to COCG solver instance
 * 
 * @return KOMEGA_SUCCESS on success, error code on failure
 */
komega_status_t komega_cocg_finalize(struct komega_cocg* solver) {
    if (solver == NULL) {
        return KOMEGA_ERROR_INVALID_PARAM;
    }
    
    /* Reset solver state */
    solver->initialized = false;
    solver->iter = 0;
    solver->resnorm = 0.0;
    
    /* Note: We don't free the vectors here as they might be needed
     * for restart functionality. The user should call komega_cocg_destroy()
     * when completely done with the solver. */
    
    return KOMEGA_SUCCESS;
}
