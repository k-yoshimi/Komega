/**
 * @file komega_vals_r.h
 * @brief Real-valued storage header for Komega C library
 */

#ifndef KOMEGA_VALS_R_H
#define KOMEGA_VALS_R_H

#include "komega.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Real values structure */
struct komega_vals_r {
    int nz;             /* Number of frequencies */
    int itermax;        /* Maximum iterations */
    double* z;          /* Frequency values */
    double* pi;         /* Current π values */
    double* pi_old;     /* Previous π values */
    double* pi_save;    /* Saved π values for restart */
    double z_seed;      /* Seed frequency */
    double rho;         /* ρ value */
    double alpha;       /* α value */
    double alpha_old;   /* Previous α value */
    double beta;        /* β value */
    bool initialized;   /* Initialization flag */
};

/* Function declarations */
komega_vals_r_t* komega_vals_r_create(void);
void komega_vals_r_destroy(komega_vals_r_t* vals);
komega_status_t komega_vals_r_initialize(komega_vals_r_t* vals, 
                                        const double* z, int nz, int itermax);

/* Value manipulation functions */
void komega_vals_r_set_z_seed(komega_vals_r_t* vals, int iz);
void komega_vals_r_scale_pi_values(komega_vals_r_t* vals, double scale_factor);
void komega_vals_r_update_pi_values(komega_vals_r_t* vals);
void komega_vals_r_save_pi_values(komega_vals_r_t* vals, int iter);
void komega_vals_r_restore_pi_values(komega_vals_r_t* vals, int iter);

/* Accessor functions */
bool komega_vals_r_is_initialized(const komega_vals_r_t* vals);
int komega_vals_r_get_nz(const komega_vals_r_t* vals);
double komega_vals_r_get_z_seed(const komega_vals_r_t* vals);
double komega_vals_r_get_rho(const komega_vals_r_t* vals);
void komega_vals_r_set_rho(komega_vals_r_t* vals, double rho);
double komega_vals_r_get_alpha(const komega_vals_r_t* vals);
void komega_vals_r_set_alpha(komega_vals_r_t* vals, double alpha);
double komega_vals_r_get_alpha_old(const komega_vals_r_t* vals);
void komega_vals_r_set_alpha_old(komega_vals_r_t* vals, double alpha_old);
double komega_vals_r_get_beta(const komega_vals_r_t* vals);
void komega_vals_r_set_beta(komega_vals_r_t* vals, double beta);
const double* komega_vals_r_get_pi_values(const komega_vals_r_t* vals);
const double* komega_vals_r_get_pi_old_values(const komega_vals_r_t* vals);

/* Global instance functions */
komega_vals_r_t* komega_get_global_vals_r(void);
void komega_cleanup_global_vals_r(void);

#ifdef __cplusplus
}
#endif

#endif /* KOMEGA_VALS_R_H */
