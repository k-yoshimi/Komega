/**
 * @file komega_vals_r.c
 * @brief Real-valued storage module for Komega C library
 * 
 * This module provides real-valued storage for CG algorithms.
 * It corresponds to the komega_vals_r Fortran module.
 */

#include "komega_vals_r.h"
#include <stdlib.h>
#include <string.h>

/* Global real values instance */
static komega_vals_r_t* global_vals_r = NULL;

komega_vals_r_t* komega_vals_r_create(void) {
    komega_vals_r_t* vals = malloc(sizeof(komega_vals_r_t));
    if (vals == NULL) {
        return NULL;
    }
    
    /* Initialize with default values */
    vals->nz = 0;
    vals->itermax = 0;
    vals->z = NULL;
    vals->pi = NULL;
    vals->pi_old = NULL;
    vals->pi_save = NULL;
    vals->z_seed = 0.0;
    vals->rho = 1.0;
    vals->alpha = 1.0;
    vals->alpha_old = 1.0;
    vals->beta = 0.0;
    vals->initialized = false;
    
    return vals;
}

void komega_vals_r_destroy(komega_vals_r_t* vals) {
    if (vals != NULL) {
        if (vals->z != NULL) {
            free(vals->z);
        }
        if (vals->pi != NULL) {
            free(vals->pi);
        }
        if (vals->pi_old != NULL) {
            free(vals->pi_old);
        }
        if (vals->pi_save != NULL) {
            free(vals->pi_save);
        }
        free(vals);
    }
}

komega_status_t komega_vals_r_initialize(komega_vals_r_t* vals, 
                                        const double* z, int nz, int itermax) {
    if (vals == NULL || z == NULL || nz <= 0 || itermax <= 0) {
        return KOMEGA_ERROR_INVALID_PARAM;
    }
    
    /* Clean up existing arrays */
    if (vals->z != NULL) {
        free(vals->z);
    }
    if (vals->pi != NULL) {
        free(vals->pi);
    }
    if (vals->pi_old != NULL) {
        free(vals->pi_old);
    }
    if (vals->pi_save != NULL) {
        free(vals->pi_save);
    }
    
    /* Allocate new arrays */
    vals->nz = nz;
    vals->itermax = itermax;
    
    vals->z = malloc(nz * sizeof(double));
    vals->pi = malloc(nz * sizeof(double));
    vals->pi_old = malloc(nz * sizeof(double));
    vals->pi_save = malloc(nz * itermax * sizeof(double));
    
    if (vals->z == NULL || vals->pi == NULL || 
        vals->pi_old == NULL || vals->pi_save == NULL) {
        komega_vals_r_destroy(vals);
        return KOMEGA_ERROR_MEMORY_ALLOCATION;
    }
    
    /* Copy input data */
    memcpy(vals->z, z, nz * sizeof(double));
    
    /* Initialize π arrays */
    for (int i = 0; i < nz; i++) {
        vals->pi[i] = 1.0;
        vals->pi_old[i] = 1.0;
    }
    
    /* Initialize π_save array */
    for (int i = 0; i < nz * itermax; i++) {
        vals->pi_save[i] = 0.0;
    }
    
    /* Initialize scalar values */
    vals->z_seed = z[0];
    vals->rho = 1.0;
    vals->alpha = 1.0;
    vals->alpha_old = 1.0;
    vals->beta = 0.0;
    vals->initialized = true;
    
    return KOMEGA_SUCCESS;
}

void komega_vals_r_set_z_seed(komega_vals_r_t* vals, int iz) {
    if (vals != NULL && vals->z != NULL && iz >= 0 && iz < vals->nz) {
        vals->z_seed = vals->z[iz];
    }
}

void komega_vals_r_scale_pi_values(komega_vals_r_t* vals, double scale_factor) {
    if (vals != NULL && vals->pi != NULL) {
        for (int i = 0; i < vals->nz; i++) {
            vals->pi[i] *= scale_factor;
        }
    }
}

void komega_vals_r_update_pi_values(komega_vals_r_t* vals) {
    if (vals != NULL && vals->pi != NULL && vals->pi_old != NULL) {
        for (int i = 0; i < vals->nz; i++) {
            vals->pi_old[i] = vals->pi[i];
        }
    }
}

void komega_vals_r_save_pi_values(komega_vals_r_t* vals, int iter) {
    if (vals != NULL && vals->pi != NULL && vals->pi_save != NULL && 
        iter >= 0 && iter < vals->itermax) {
        for (int i = 0; i < vals->nz; i++) {
            vals->pi_save[i * vals->itermax + iter] = vals->pi[i];
        }
    }
}

void komega_vals_r_restore_pi_values(komega_vals_r_t* vals, int iter) {
    if (vals != NULL && vals->pi != NULL && vals->pi_save != NULL && 
        iter >= 0 && iter < vals->itermax) {
        for (int i = 0; i < vals->nz; i++) {
            vals->pi[i] = vals->pi_save[i * vals->itermax + iter];
        }
    }
}

/* Accessor functions */
bool komega_vals_r_is_initialized(const komega_vals_r_t* vals) {
    return (vals != NULL) ? vals->initialized : false;
}

int komega_vals_r_get_nz(const komega_vals_r_t* vals) {
    return (vals != NULL) ? vals->nz : 0;
}

double komega_vals_r_get_z_seed(const komega_vals_r_t* vals) {
    return (vals != NULL) ? vals->z_seed : 0.0;
}

double komega_vals_r_get_rho(const komega_vals_r_t* vals) {
    return (vals != NULL) ? vals->rho : 0.0;
}

void komega_vals_r_set_rho(komega_vals_r_t* vals, double rho) {
    if (vals != NULL) {
        vals->rho = rho;
    }
}

double komega_vals_r_get_alpha(const komega_vals_r_t* vals) {
    return (vals != NULL) ? vals->alpha : 0.0;
}

void komega_vals_r_set_alpha(komega_vals_r_t* vals, double alpha) {
    if (vals != NULL) {
        vals->alpha = alpha;
    }
}

double komega_vals_r_get_alpha_old(const komega_vals_r_t* vals) {
    return (vals != NULL) ? vals->alpha_old : 0.0;
}

void komega_vals_r_set_alpha_old(komega_vals_r_t* vals, double alpha_old) {
    if (vals != NULL) {
        vals->alpha_old = alpha_old;
    }
}

double komega_vals_r_get_beta(const komega_vals_r_t* vals) {
    return (vals != NULL) ? vals->beta : 0.0;
}

void komega_vals_r_set_beta(komega_vals_r_t* vals, double beta) {
    if (vals != NULL) {
        vals->beta = beta;
    }
}

const double* komega_vals_r_get_pi_values(const komega_vals_r_t* vals) {
    return (vals != NULL) ? vals->pi : NULL;
}

const double* komega_vals_r_get_pi_old_values(const komega_vals_r_t* vals) {
    return (vals != NULL) ? vals->pi_old : NULL;
}

komega_vals_r_t* komega_get_global_vals_r(void) {
    if (global_vals_r == NULL) {
        global_vals_r = komega_vals_r_create();
    }
    return global_vals_r;
}

void komega_cleanup_global_vals_r(void) {
    if (global_vals_r != NULL) {
        komega_vals_r_destroy(global_vals_r);
        global_vals_r = NULL;
    }
}
