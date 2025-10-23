/**
 * @file komega_fortran_wrapper.h
 * @brief Fortran library wrapper header for C language
 * 
 * This header provides C language wrappers for the existing Fortran Komega library.
 * It allows C programs to call Fortran functions directly.
 */

#ifndef KOMEGA_FORTRAN_WRAPPER_H
#define KOMEGA_FORTRAN_WRAPPER_H

#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif

/* BiCG Solver Functions */
void komega_bicg_init(int *ndim, int *nl, int *nz, double complex *x,
                      double complex *z, int *itermax, double *threshold, int *comm);
void komega_bicg_restart(int *ndim, int *nl, int *nz, double complex *x,
                         double complex *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double complex *v2,
                         double complex *v12, double complex *v4, double complex *v14,
                         double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save, int *comm);
void komega_bicg_update(double complex *v12, double complex *v2,
                        double complex *v14, double complex *v4,
                        double complex *x, double complex *r_l, int *status);
void komega_bicg_getcoef(double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save);
void komega_bicg_getvec(double complex *r_old, double complex *r_tilde_old);
void komega_bicg_getresidual(double *res);
void komega_bicg_finalize(void);

/* COCG Solver Functions */
void komega_cocg_init(int *ndim, int *nl, int *nz, double complex *x,
                      double complex *z, int *itermax, double *threshold, int *comm);
void komega_cocg_restart(int *ndim, int *nl, int *nz, double complex *x,
                         double complex *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double complex *v2,
                         double complex *v12, 
                         double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save, int *comm);
void komega_cocg_update(double complex *v12, double complex *v2,
                        double complex *x, double complex *r_l, int *status);
void komega_cocg_getcoef(double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save);
void komega_cocg_getvec(double complex *r_old);
void komega_cocg_getresidual(double *res);
void komega_cocg_finalize(void);

/* CG Complex Solver Functions */
void komega_cg_c_init(int *ndim, int *nl, int *nz, double complex *x,
                      double *z, int *itermax, double *threshold, int *comm);
void komega_cg_c_restart(int *ndim, int *nl, int *nz, double complex *x,
                         double *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double complex *v2,
                         double complex *v12, 
                         double *alpha_save, double *beta_save,
                         double *z_seed, double complex *r_l_save, int *comm);
void komega_cg_c_update(double complex *v12, double complex *v2,
                        double complex *x, double complex *r_l, int *status);
void komega_cg_c_getcoef(double *alpha_save, double *beta_save,
                         double *z_seed, double complex *r_l_save);
void komega_cg_c_getvec(double complex *r_old);
void komega_cg_c_getresidual(double *res);
void komega_cg_c_finalize(void);

/* CG Real Solver Functions */
void komega_cg_r_init(int *ndim, int *nl, int *nz, double *x,
                      double *z, int *itermax, double *threshold, int *comm);
void komega_cg_r_restart(int *ndim, int *nl, int *nz, double *x,
                         double *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double *v2,
                         double *v12, 
                         double *alpha_save, double *beta_save,
                         double *z_seed, double *r_l_save, int *comm);
void komega_cg_r_update(double *v12, double *v2,
                        double *x, double *r_l, int *status);
void komega_cg_r_getcoef(double *alpha_save, double *beta_save,
                         double *z_seed, double *r_l_save);
void komega_cg_r_getvec(double *r_old);
void komega_cg_r_getresidual(double *res);
void komega_cg_r_finalize(void);

#ifdef __cplusplus
}
#endif

#endif /* KOMEGA_FORTRAN_WRAPPER_H */
