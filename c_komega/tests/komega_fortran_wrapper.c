/**
 * @file komega_fortran_wrapper.c
 * @brief Fortran library wrapper for C language
 * 
 * This module provides C language wrappers for the existing Fortran Komega library.
 * It allows C programs to call Fortran functions directly.
 */

#include "komega_fortran_wrapper.h"
#include <stdlib.h>
#include <string.h>

/* Fortran function declarations (from src/komega.h) */
extern void komega_bicg_init_(int *ndim, int *nl, int *nz, double complex *x,
                              double complex *z, int *itermax, double *threshold, int *comm);
extern void komega_cocg_init_(int *ndim, int *nl, int *nz, double complex *x,
                              double complex *z, int *itermax, double *threshold, int *comm);
extern void komega_cg_c_init_(int *ndim, int *nl, int *nz, double complex *x,
                              double *z, int *itermax, double *threshold, int *comm);
extern void komega_cg_r_init_(int *ndim, int *nl, int *nz, double *x,
                              double *z, int *itermax, double *threshold, int *comm);

extern void komega_bicg_restart_(int *ndim, int *nl, int *nz, double complex *x,
                                 double complex *z, int *itermax, double *threshold,
                                 int *status, int *iter_old, double complex *v2,
                                 double complex *v12, double complex *v4, double complex *v14,
                                 double complex *alpha_save, double complex *beta_save,
                                 double complex *z_seed, double complex *r_l_save, int *comm);
extern void komega_cocg_restart_(int *ndim, int *nl, int *nz, double complex *x,
                                 double complex *z, int *itermax, double *threshold,
                                 int *status, int *iter_old, double complex *v2,
                                 double complex *v12, 
                                 double complex *alpha_save, double complex *beta_save,
                                 double complex *z_seed, double complex *r_l_save, int *comm);
extern void komega_cg_c_restart_(int *ndim, int *nl, int *nz, double complex *x,
                                 double *z, int *itermax, double *threshold,
                                 int *status, int *iter_old, double complex *v2,
                                 double complex *v12, 
                                 double *alpha_save, double *beta_save,
                                 double *z_seed, double complex *r_l_save, int *comm);
extern void komega_cg_r_restart_(int *ndim, int *nl, int *nz, double *x,
                                 double *z, int *itermax, double *threshold,
                                 int *status, int *iter_old, double *v2,
                                 double *v12, 
                                 double *alpha_save, double *beta_save,
                                 double *z_seed, double *r_l_save, int *comm);

extern void komega_bicg_update_(double complex *v12, double complex *v2,
                                double complex *v14, double complex *v4,
                                double complex *x, double complex *r_l, int *status);
extern void komega_cocg_update_(double complex *v12, double complex *v2,
                                double complex *x, double complex *r_l, int *status);
extern void komega_cg_c_update_(double complex *v12, double complex *v2,
                                double complex *x, double complex *r_l, int *status);
extern void komega_cg_r_update_(double *v12, double *v2,
                                double *x, double *r_l, int *status);

extern void komega_bicg_getcoef_(double complex *alpha_save, double complex *beta_save,
                                 double complex *z_seed, double complex *r_l_save);
extern void komega_cocg_getcoef_(double complex *alpha_save, double complex *beta_save,
                                 double complex *z_seed, double complex *r_l_save);
extern void komega_cg_c_getcoef_(double *alpha_save, double *beta_save,
                                 double *z_seed, double complex *r_l_save);
extern void komega_cg_r_getcoef_(double *alpha_save, double *beta_save,
                                 double *z_seed, double *r_l_save);

extern void komega_bicg_getvec_(double complex *r_old, double complex *r_tilde_old);
extern void komega_cocg_getvec_(double complex *r_old);
extern void komega_cg_c_getvec_(double complex *r_old);
extern void komega_cg_r_getvec_(double *r_old);

extern void komega_bicg_getresidual_(double *res);
extern void komega_cocg_getresidual_(double *res);
extern void komega_cg_c_getresidual_(double *res);
extern void komega_cg_r_getresidual_(double *res);

extern void komega_bicg_finalize_(void);
extern void komega_cocg_finalize_(void);
extern void komega_cg_r_finalize_(void);
extern void komega_cg_c_finalize_(void);

/* C wrapper functions for BiCG solver */
void komega_bicg_init(int *ndim, int *nl, int *nz, double complex *x,
                      double complex *z, int *itermax, double *threshold, int *comm) {
    komega_bicg_init_(ndim, nl, nz, x, z, itermax, threshold, comm);
}

void komega_bicg_restart(int *ndim, int *nl, int *nz, double complex *x,
                         double complex *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double complex *v2,
                         double complex *v12, double complex *v4, double complex *v14,
                         double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save, int *comm) {
    komega_bicg_restart_(ndim, nl, nz, x, z, itermax, threshold, status, iter_old, v2,
                        v12, v4, v14, alpha_save, beta_save, z_seed, r_l_save, comm);
}

void komega_bicg_update(double complex *v12, double complex *v2,
                        double complex *v14, double complex *v4,
                        double complex *x, double complex *r_l, int *status) {
    komega_bicg_update_(v12, v2, v14, v4, x, r_l, status);
}

void komega_bicg_getcoef(double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save) {
    komega_bicg_getcoef_(alpha_save, beta_save, z_seed, r_l_save);
}

void komega_bicg_getvec(double complex *r_old, double complex *r_tilde_old) {
    komega_bicg_getvec_(r_old, r_tilde_old);
}

void komega_bicg_getresidual(double *res) {
    komega_bicg_getresidual_(res);
}

void komega_bicg_finalize(void) {
    komega_bicg_finalize_();
}

/* C wrapper functions for COCG solver */
void komega_cocg_init(int *ndim, int *nl, int *nz, double complex *x,
                       double complex *z, int *itermax, double *threshold, int *comm) {
    komega_cocg_init_(ndim, nl, nz, x, z, itermax, threshold, comm);
}

void komega_cocg_restart(int *ndim, int *nl, int *nz, double complex *x,
                         double complex *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double complex *v2,
                         double complex *v12, 
                         double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save, int *comm) {
    komega_cocg_restart_(ndim, nl, nz, x, z, itermax, threshold, status, iter_old, v2,
                        v12, alpha_save, beta_save, z_seed, r_l_save, comm);
}

void komega_cocg_update(double complex *v12, double complex *v2,
                        double complex *x, double complex *r_l, int *status) {
    komega_cocg_update_(v12, v2, x, r_l, status);
}

void komega_cocg_getcoef(double complex *alpha_save, double complex *beta_save,
                         double complex *z_seed, double complex *r_l_save) {
    komega_cocg_getcoef_(alpha_save, beta_save, z_seed, r_l_save);
}

void komega_cocg_getvec(double complex *r_old) {
    komega_cocg_getvec_(r_old);
}

void komega_cocg_getresidual(double *res) {
    komega_cocg_getresidual_(res);
}

void komega_cocg_finalize(void) {
    komega_cocg_finalize_();
}

/* C wrapper functions for CG Complex solver */
void komega_cg_c_init(int *ndim, int *nl, int *nz, double complex *x,
                       double *z, int *itermax, double *threshold, int *comm) {
    komega_cg_c_init_(ndim, nl, nz, x, z, itermax, threshold, comm);
}

void komega_cg_c_restart(int *ndim, int *nl, int *nz, double complex *x,
                         double *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double complex *v2,
                         double complex *v12, 
                         double *alpha_save, double *beta_save,
                         double *z_seed, double complex *r_l_save, int *comm) {
    komega_cg_c_restart_(ndim, nl, nz, x, z, itermax, threshold, status, iter_old, v2,
                        v12, alpha_save, beta_save, z_seed, r_l_save, comm);
}

void komega_cg_c_update(double complex *v12, double complex *v2,
                        double complex *x, double complex *r_l, int *status) {
    komega_cg_c_update_(v12, v2, x, r_l, status);
}

void komega_cg_c_getcoef(double *alpha_save, double *beta_save,
                         double *z_seed, double complex *r_l_save) {
    komega_cg_c_getcoef_(alpha_save, beta_save, z_seed, r_l_save);
}

void komega_cg_c_getvec(double complex *r_old) {
    komega_cg_c_getvec_(r_old);
}

void komega_cg_c_getresidual(double *res) {
    komega_cg_c_getresidual_(res);
}

void komega_cg_c_finalize(void) {
    komega_cg_c_finalize_();
}

/* C wrapper functions for CG Real solver */
void komega_cg_r_init(int *ndim, int *nl, int *nz, double *x,
                       double *z, int *itermax, double *threshold, int *comm) {
    komega_cg_r_init_(ndim, nl, nz, x, z, itermax, threshold, comm);
}

void komega_cg_r_restart(int *ndim, int *nl, int *nz, double *x,
                         double *z, int *itermax, double *threshold,
                         int *status, int *iter_old, double *v2,
                         double *v12, 
                         double *alpha_save, double *beta_save,
                         double *z_seed, double *r_l_save, int *comm) {
    komega_cg_r_restart_(ndim, nl, nz, x, z, itermax, threshold, status, iter_old, v2,
                        v12, alpha_save, beta_save, z_seed, r_l_save, comm);
}

void komega_cg_r_update(double *v12, double *v2,
                       double *x, double *r_l, int *status) {
    komega_cg_r_update_(v12, v2, x, r_l, status);
}

void komega_cg_r_getcoef(double *alpha_save, double *beta_save,
                         double *z_seed, double *r_l_save) {
    komega_cg_r_getcoef_(alpha_save, beta_save, z_seed, r_l_save);
}

void komega_cg_r_getvec(double *r_old) {
    komega_cg_r_getvec_(r_old);
}

void komega_cg_r_getresidual(double *res) {
    komega_cg_r_getresidual_(res);
}

void komega_cg_r_finalize(void) {
    komega_cg_r_finalize_();
}
