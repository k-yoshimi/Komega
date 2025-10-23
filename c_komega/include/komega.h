/**
 * @file komega.h
 * @brief Main header file for the Komega C library
 * 
 * This file provides the main interface for the Komega C library,
 * which is a C port of the Fortran Komega library for solving
 * linear systems using iterative methods.
 * 
 * The Komega library provides iterative solvers for linear systems,
 * including BiCG, CG, and COCG algorithms. It is designed to be
 * compatible with the original Fortran implementation while providing
 * a modern C interface.
 * 
 * @author Mitsuaki Kawamura (Original Fortran)
 * @author C Port for verification and testing purposes
 * @date 2025
 * @version 1.0.0
 * 
 * @see https://github.com/k-yoshimi/Komega
 * @see src/komega_parameter.h
 * @see src/komega_math.h
 * @see src/komega_vals_r.h
 */

#ifndef KOMEGA_H
#define KOMEGA_H

#ifdef __cplusplus
extern "C" {
#endif

#include <complex.h>
#include <stdbool.h>
#include <stddef.h>

/**
 * @defgroup Version Version Information
 * @{
 */

/** @brief Major version number */
#define KOMEGA_VERSION_MAJOR 1

/** @brief Minor version number */
#define KOMEGA_VERSION_MINOR 0

/** @brief Patch version number */
#define KOMEGA_VERSION_PATCH 0

/** @brief Version string */
#define KOMEGA_VERSION_STRING "1.0.0"

/** @} */

/**
 * @defgroup Constants Mathematical Constants
 * @{
 */

/** @brief Pi constant */
#define KOMEGA_PI 3.14159265358979323846

/** @brief Machine epsilon */
#define KOMEGA_EPSILON 1e-15

/** @brief Almost zero threshold */
#define KOMEGA_ALMOST_ZERO 1e-12

/** @} */

/**
 * @defgroup SolverTypes Solver Types
 * @{
 */

/** @brief Bi-Conjugate Gradient solver */
typedef enum {
    KOMEGA_SOLVER_BICG = 0,    /**< BiCG solver for complex systems */
    KOMEGA_SOLVER_CG_R,        /**< CG solver for real systems */
    KOMEGA_SOLVER_CG_C,        /**< CG solver for complex systems */
    KOMEGA_SOLVER_COCG         /**< COCG solver for complex systems */
} komega_solver_type_t;

/** @} */

/**
 * @defgroup StatusCodes Status Codes
 * @{
 */

/** @brief Status codes for Komega operations */
typedef enum {
    KOMEGA_SUCCESS = 0,                    /**< Operation successful */
    KOMEGA_ERROR_INVALID_PARAM,           /**< Invalid parameter */
    KOMEGA_ERROR_MEMORY_ALLOCATION,       /**< Memory allocation failed */
    KOMEGA_ERROR_SOLVER_NOT_INITIALIZED,  /**< Solver not initialized */
    KOMEGA_ERROR_CONVERGENCE_FAILED,      /**< Convergence failed */
    KOMEGA_ERROR_UNKNOWN_SOLVER_TYPE      /**< Unknown solver type */
} komega_status_t;

/** @} */

/* Complex number type */
typedef double complex komega_complex_t;

/* Forward declarations */
typedef struct komega_parameter komega_parameter_t;
typedef struct komega_math komega_math_t;
typedef struct komega_vals_r komega_vals_r_t;
typedef struct komega_vals_c komega_vals_c_t;
typedef struct komega_vecs_r komega_vecs_r_t;
typedef struct komega_vecs_c komega_vecs_c_t;
typedef struct komega_solver komega_solver_t;

/* Solver type declarations */
typedef struct komega_bicg komega_bicg_t;
typedef struct komega_cg_r komega_cg_r_t;
typedef struct komega_cg_c komega_cg_c_t;
typedef struct komega_cocg komega_cocg_t;

/* Parameter module functions */
komega_parameter_t* komega_parameter_create(void);
void komega_parameter_destroy(komega_parameter_t* params);
komega_status_t komega_parameter_initialize(komega_parameter_t* params, 
                                           int ndim, int nl, int nz, 
                                           int itermax, double threshold);

/* Math module functions */
komega_math_t* komega_math_create(void);
void komega_math_destroy(komega_math_t* math);
double komega_math_ddot(const double* x, const double* y, int n);
komega_complex_t komega_math_zdotc(const komega_complex_t* x, 
                                   const komega_complex_t* y, int n);
komega_complex_t komega_math_zdotu(const komega_complex_t* x, 
                                   const komega_complex_t* y, int n);
void komega_math_dscal(double alpha, double* x, int n);
void komega_math_zscal(komega_complex_t alpha, komega_complex_t* x, int n);
void komega_math_daxpy(double alpha, const double* x, double* y, int n);
void komega_math_zaxpy(komega_complex_t alpha, const komega_complex_t* x, 
                       komega_complex_t* y, int n);
void komega_math_dcopy(const double* x, double* y, int n);
void komega_math_zcopy(const komega_complex_t* x, komega_complex_t* y, int n);
double komega_math_dabsmax(const double* x, int n);
double komega_math_zabsmax(const komega_complex_t* x, int n);

/* Value storage functions */
komega_vals_r_t* komega_vals_r_create(void);
void komega_vals_r_destroy(komega_vals_r_t* vals);
komega_status_t komega_vals_r_initialize(komega_vals_r_t* vals, 
                                        const double* z, int nz, int itermax);

komega_vals_c_t* komega_vals_c_create(void);
void komega_vals_c_destroy(komega_vals_c_t* vals);
komega_status_t komega_vals_c_initialize(komega_vals_c_t* vals, 
                                        const komega_complex_t* z, int nz, int itermax);

/* Vector storage functions */
komega_vecs_r_t* komega_vecs_r_create(void);
void komega_vecs_r_destroy(komega_vecs_r_t* vecs);
komega_status_t komega_vecs_r_initialize(komega_vecs_r_t* vecs, 
                                         int ndim, int nl, int nz, int itermax);

komega_vecs_c_t* komega_vecs_c_create(void);
void komega_vecs_c_destroy(komega_vecs_c_t* vecs);
komega_status_t komega_vecs_c_initialize(komega_vecs_c_t* vecs, 
                                         int ndim, int nl, int nz, int itermax);

/* Solver functions */
komega_solver_t* komega_solver_create(komega_solver_type_t type);
void komega_solver_destroy(komega_solver_t* solver);
komega_status_t komega_solver_init(komega_solver_t* solver, 
                                   int ndim, int nl, int nz,
                                   const void* z_frequencies, 
                                   int itermax, double threshold);
komega_status_t komega_solver_update(komega_solver_t* solver, 
                                     const void* v12, const void* v2, 
                                     void* x, void* r_l, int* status);
komega_status_t komega_solver_get_coefficients(komega_solver_t* solver, 
                                               void* alpha_save, void* beta_save,
                                               void* z_seed, void* r_l_save);
komega_status_t komega_solver_get_vectors(komega_solver_t* solver, void* vectors);
komega_status_t komega_solver_get_residual(komega_solver_t* solver, double* residual);
komega_status_t komega_solver_finalize(komega_solver_t* solver);

/* Utility functions */
const char* komega_get_version_string(void);
const char* komega_get_solver_name(komega_solver_type_t type);
const char* komega_get_status_message(komega_status_t status);
bool komega_solver_supports_restart(komega_solver_type_t type);
bool komega_solver_supports_mpi(komega_solver_type_t type);
const char* komega_solver_get_data_type(komega_solver_type_t type);

/* Global instances (for compatibility with Fortran/Python versions) */
komega_parameter_t* komega_get_global_parameter(void);
komega_math_t* komega_get_global_math(void);
komega_vals_r_t* komega_get_global_vals_r(void);
komega_vals_c_t* komega_get_global_vals_c(void);
komega_vecs_r_t* komega_get_global_vecs_r(void);
komega_vecs_c_t* komega_get_global_vecs_c(void);

/* Solver creation functions */
komega_bicg_t* komega_bicg_create(void);
komega_cg_r_t* komega_cg_r_create(void);
komega_cg_c_t* komega_cg_c_create(void);
komega_cocg_t* komega_cocg_create(void);

/* Solver destruction functions */
void komega_bicg_destroy(komega_bicg_t* solver);
void komega_cg_r_destroy(komega_cg_r_t* solver);
void komega_cg_c_destroy(komega_cg_c_t* solver);
void komega_cocg_destroy(komega_cocg_t* solver);

#ifdef __cplusplus
}
#endif

#endif /* KOMEGA_H */
