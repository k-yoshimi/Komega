"""
Komega Python Library - A library for solving linear systems in materials science

This is a Python port of the Komega Fortran library, providing the same
functionality for solving linear systems with multiple frequency points
using various iterative methods.

Copyright (C) 2016 Mitsuaki Kawamura
Python port created for verification and testing purposes.
"""

from .komega import (
    KomegaSolver,
    create_solver,
    get_available_solvers,
    get_bicg_solver,
    get_cg_c_solver,
    get_cg_r_solver,
    get_cocg_solver,
    get_solver_info,
)
from .komega_bicg import (
    KomegaBiCG,
    get_global_bicg,
    komega_BICG_finalize,
    komega_BICG_getcoef,
    komega_BICG_getresidual,
    komega_BICG_getvec,
    komega_BICG_init,
    komega_BICG_update,
)
from .komega_cg_c import (
    KomegaCGC,
    get_global_cg_c,
    komega_CG_C_finalize,
    komega_CG_C_getcoef,
    komega_CG_C_getresidual,
    komega_CG_C_getvec,
    komega_CG_C_init,
    komega_CG_C_update,
)
from .komega_cg_r import (
    KomegaCGR,
    get_global_cg_r,
    komega_CG_R_finalize,
    komega_CG_R_getcoef,
    komega_CG_R_getresidual,
    komega_CG_R_getvec,
    komega_CG_R_init,
    komega_CG_R_update,
)
from .komega_cocg import (
    KomegaCOCG,
    get_global_cocg,
    komega_COCG_finalize,
    komega_COCG_getcoef,
    komega_COCG_getresidual,
    komega_COCG_getvec,
    komega_COCG_init,
    komega_COCG_update,
)
from .komega_math import (
    KomegaMath,
    dabsmax,
    daxpy,
    dcopy,
    ddot,
    ddotMPI,
    dscal,
    get_global_math,
    zaxpy,
    zcopy,
    zdotc,
    zdotcMPI,
    zdotu,
    zdotuMPI,
    zscal,
)
from .komega_parameter import (
    KomegaParameter,
    get_global_params,
    initialize_global_params,
    reset_global_params,
)
from .komega_vals_c import (
    KomegaValsC,
    cleanup_vals_c,
    get_global_vals_c,
    initialize_vals_c,
)
from .komega_vals_r import (
    KomegaValsR,
    cleanup_vals_r,
    get_global_vals_r,
    initialize_vals_r,
)
from .komega_vecs_c import (
    KomegaVecsC,
    cleanup_vecs_c,
    get_global_vecs_c,
    initialize_vecs_c,
)
from .komega_vecs_r import (
    KomegaVecsR,
    cleanup_vecs_r,
    get_global_vecs_r,
    initialize_vecs_r,
)

# Version information
__version__ = "2.0.0"
__author__ = "Mitsuaki Kawamura"
__email__ = "kawamura@issp.u-tokyo.ac.jp"
__description__ = "Komega - A library for solving linear systems in materials science"

# Main classes and functions
__all__ = [
    # Main interface
    "KomegaSolver",
    "create_solver",
    "get_available_solvers",
    "get_solver_info",
    # Global solvers
    "get_bicg_solver",
    "get_cg_r_solver",
    "get_cg_c_solver",
    "get_cocg_solver",
    # Parameter management
    "KomegaParameter",
    "get_global_params",
    "initialize_global_params",
    "reset_global_params",
    # Math operations
    "KomegaMath",
    "get_global_math",
    "ddot",
    "zdotc",
    "zdotu",
    "dscal",
    "zscal",
    "dcopy",
    "zcopy",
    "daxpy",
    "zaxpy",
    "ddotMPI",
    "zdotcMPI",
    "zdotuMPI",
    "dabsmax",
    # Value storage
    "KomegaValsR",
    "get_global_vals_r",
    "initialize_vals_r",
    "cleanup_vals_r",
    "KomegaValsC",
    "get_global_vals_c",
    "initialize_vals_c",
    "cleanup_vals_c",
    # Vector storage
    "KomegaVecsR",
    "get_global_vecs_r",
    "initialize_vecs_r",
    "cleanup_vecs_r",
    "KomegaVecsC",
    "get_global_vecs_c",
    "initialize_vecs_c",
    "cleanup_vecs_c",
    # Solvers
    "KomegaBiCG",
    "get_global_bicg",
    "komega_BICG_init",
    "komega_BICG_update",
    "komega_BICG_getcoef",
    "komega_BICG_getvec",
    "komega_BICG_getresidual",
    "komega_BICG_finalize",
    "KomegaCGR",
    "get_global_cg_r",
    "komega_CG_R_init",
    "komega_CG_R_update",
    "komega_CG_R_getcoef",
    "komega_CG_R_getvec",
    "komega_CG_R_getresidual",
    "komega_CG_R_finalize",
    "KomegaCGC",
    "get_global_cg_c",
    "komega_CG_C_init",
    "komega_CG_C_update",
    "komega_CG_C_getcoef",
    "komega_CG_C_getvec",
    "komega_CG_C_getresidual",
    "komega_CG_C_finalize",
    "KomegaCOCG",
    "get_global_cocg",
    "komega_COCG_init",
    "komega_COCG_update",
    "komega_COCG_getcoef",
    "komega_COCG_getvec",
    "komega_COCG_getresidual",
    "komega_COCG_finalize",
]
