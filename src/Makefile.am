#
# Komega library
#
MAINTAINERCLEANFILES = Makefile.in

include_HEADERS = komega.h
# Remove nodist_include_HEADERS to avoid dependency issues
# nodist_include_HEADERS = komega_bicg.$(FC_MODEXT) komega_cg_c.$(FC_MODEXT) \
#						 komega_cg_r.$(FC_MODEXT) komega_cocg.$(FC_MODEXT) \
#						 komega_shifted_qmr_sym.$(FC_MODEXT) \
#						 komega_shifted_qmr_sym_b.$(FC_MODEXT)

lib_LTLIBRARIES = libkomega.la

libkomega_la_SOURCES = komega_bicg.F90 komega_cg_r.F90 komega_cg_c.F90 komega_cocg.F90 \
                       komega_shifted_qmr_sym.F90 komega_shifted_qmr_sym_b.F90 komega_math.F90 \
                       komega_vals.F90

# Ensure Fortran modules are built before other targets
BUILT_SOURCES = komega_vals.$(FC_MODEXT) komega_math.$(FC_MODEXT) \
               komega_bicg.$(FC_MODEXT) komega_cg_c.$(FC_MODEXT) \
               komega_cg_r.$(FC_MODEXT) komega_cocg.$(FC_MODEXT) \
               komega_shifted_qmr_sym.$(FC_MODEXT) \
               komega_shifted_qmr_sym_b.$(FC_MODEXT)

# Clean up modules.
clean-local: clean-modules
clean-modules:
	test -z "$(FC_MODEXT)" || $(RM) *.$(FC_MODEXT)

# Fix dependency order - build komega_vals first, then komega_math, then others
komega_vals.lo : komega_vals.F90
komega_math.lo : komega_math.F90 komega_vals.lo
komega_bicg.lo : komega_bicg.F90 komega_math.lo komega_vals.lo
komega_cg_c.lo : komega_cg_c.F90 komega_math.lo komega_vals.lo
komega_cg_r.lo : komega_cg_r.F90 komega_math.lo komega_vals.lo
komega_cocg.lo : komega_cocg.F90 komega_math.lo komega_vals.lo
komega_shifted_qmr_sym.lo : komega_shifted_qmr_sym.F90 komega_math.lo komega_vals.lo
komega_shifted_qmr_sym_b.lo : komega_shifted_qmr_sym_b.F90 komega_math.lo komega_vals.lo 

