################################################################################
# Parallel SELFE Makefile
#
# User makes environment settings for particular OS / PLATFORM / COMPILER / MPI
# below. These settings are based on the following options.
#
# Compiler name:
#   FCS: Serial compiler (for utilities)
#   FCP: Parallel compiler
#   FLD: Linker (in general same as parallel compiler)
#
# Compilation flags
#   FCSFLAGS: Flags for serial compilation
#   FCPFLAGS: Flags for parallel compilation (including all pre-processing flags)
#   FLDFLAGS: Flags for linker (e.g., -O2)
#
# Preprocessor flags:
#   DEBUG: Enable debugging code
#   ORDERED_SUM: Enable globally ordered sums & dot-products for bit reproducibility
#     of state quantities independent of number of processors (note: this can
#     significantly degrade performance);
#   INCLUDE_TIMING: Enable wallclock timing of code (note: this can have slight
#     effect on performance);
#   MPI_VERSION = 1 or 2: Version of MPI (try 2 first, if compile fails due to mpi
#     related errors then switch to version 1;
#
# Libraries (needed for parallel code). Need to re-compile these when compilers are changed.
#   MTSLIBS: Flags for linking ParMeTiS/MeTiS libaries
#   netcdf: v3.x has been tested
#   GOTM:    v3.2.5
#   ALTLIBS: (obsolete) Flags for linking alternate solver libraries (LAPACK or ITPACK,
#            these are just for testing)
################################################################################


########### Ensure that the user's Make.defs.local exists and include it

make_defs_local_exists := $(wildcard ../mk/Make.defs.local)
ifeq ($(strip $(make_defs_local_exists)),)
   $(error Make.defs.local does not exist in the /mk directory (parallel to this /src directory). See readme)
endif


-include ../mk/Make.defs.local  #optional include
include ../mk/Make.defs.config

ifeq ($(findstring GNU Fortran,$(shell $(FCP) --version)),GNU Fortran)
  $(print Here)
  MODULEOUTDIR = -J
else
  ifeq ($(ENV),DKRZ) ##please edit here, Eli
    MODULEOUTDIR = -qmoddir=
  else
    MODULEOUTDIR = -module 
  endif
endif


##########
# Check for serial compiler and revert to parallel if not provided
#########
ifeq ($(FCS),)
    FCS := ${FCP}
endif



############
# Derived compiler flags
###########

ifeq ($(ENV),DKRZ)
  IBMCPP = -WF,
else
  IBMCPP =
endif

PPFLAGS := $(PPFLAGS) $(IBMCPP)-DSELFE 

ifneq ($(USE_OPEN64),)
    ifeq ($(USE_64), yes)
      AVOID_ADV_WRITE := $(USE_OPEN64)
      PPFLAGS := $(PPFLAGS) $(IBMCPP)-DAVOID_ADV_WRITE $(IBMCPP)-DUSE_OPEN64
    endif
endif

################################################################################
# Object files; modules compile first; modules that depend on other modules compile last.
################################################################################

# rules to build objects for executable programs and libraries

#Search path
VPATH = Core/ Driver/ Hydro/ EcoSim/ ICM/ Sediment/ Sed2d/ WWMIII/ \
        /disk/ambcs03/0/yinglong/SELFE/from_Durski_Sept2009/V3.0c/ \
        TIMOR/

ifdef USE_PETSC
  include ${PETSC_DIR}/conf/variables
endif

OBJS = elfe_glbl.o elfe_msgp.o misc_modules.o selfe_driver.o

ifdef CHOOSE_TVD
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DCHOOSE_TVD
endif

ifdef USE_DEBUG
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DDEBUG 
endif

ifdef USE_TIMER
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DINCLUDE_TIMING
endif

ifdef USE_TIMOR
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_TIMOR
  OBJS := $(OBJS) flmud_input.o flmud_main.o flmud_pool.o
endif

ifdef USE_GOTM
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_GOTM
endif

ifdef USE_HA
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_HA
  OBJS := $(OBJS) harm.o
endif

ifdef USE_ECO
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_ECO
  OBJS := $(OBJS) bio_param.o eclight.o biology.o bio_init.o ecosim.o read_ecoin.o specir.o
endif

ifdef USE_ICM
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_ICM
  OBJS := $(OBJS) icm_mod.o icm_sed_mod.o icm_init.o icm_biology.o  icm.o  read_icm_input.o
endif

ifdef USE_NAPZD
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_NAPZD
  OBJS := $(OBJS) biology_napzd.o napzd_spitz.o read_napzd_input.o
endif

ifdef USE_SED
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_SED
  OBJS := $(OBJS) sed_mod.o sed_misc_subs.o \
                  read_sed_input.o sed_init.o  \
                  sed_friction.o sed_filter.o \
                  sediment.o sed_bedload.o

endif

ifdef USE_OIL
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_OIL
endif

ifdef USE_WWM
# please take care we have now ST41 (ardhuin_old) and ST42 (ardhuin_new)
  PPFLAGS := $(PPFLAGS) -DUSE_WWM -DNCDF -DST42 -DWWM_SOLVER 
OBJS := $(OBJS) wwm_datapl.o wwm_petscpool.o  wwm_petsc_seriell.o wwm_petsc_parallel.o \
       wwm_petsc_block.o wwm_petsc_controller.o wwm_aux.o wwm_aux_parall.o wwm_mjdv2.o wwm_blasaux.o wwm_sparskit.o \
       wwm_ardhuin_old.o wwm_wind.o wwm_ardhuin_new.o wwm_breaking.o wwm_friction.o wwm_cycle3.o wwm_dislin.o wwm_diclass.o \
       wwm_gridcf.o wwm_hotfile.o wwm_parall_solver.o wwm_m_constants.o wwm_m_fileio.o wwm_serv_xnl4v5.o wwm_mod_xnl4v5.o wwm_snl4_wrt.o \
       wwm_wave_setup.o wwm_initio.o wwm_netcdf.o wwm_input.o wwm_bdcons.o \
       wwm_specparam.o wwm_windinput.o wwm_dissipation.o wwm_snl3.o wwm_snl4.o wwm_snl42.o wwm_babanin.o wwm_sourceterms.o \
       wwm_specint.o wwm_nums1d.o wwm_numtheta.o wwm_numsigma.o wwm_fluctsplit.o \
       wwm_snonlin.o wwm_snonlin_local.o wwm_stress.o wwm_stresso.o wwm_stresso_local.o wwm_sbottom.o wwm_sdiss_ardh_vec.o wwm_sdiss_ardh_vec_local.o \
       wwm_sinput.o wwm_sinput_local.o wwm_sinput_ard.o wwm_sinput_ard_local.o wwm_wsigstar.o wwm_wsigstar_local.o wwm_tauhf.o wwm_airsea.o wwm_airsea_local.o \
       wwm_jafu.o wwm_nlweigt.o wwm_buildstress.o wwm_inisnonlin.o wwm_fkmean.o wwm_fkmean_local.o wwm_frcutindex.o wwm_frcutindex_local.o wwm_femeanws.o \
       wwm_femeanws_local.o wwm_sdissip.o wwm_sdissip_local.o wwm_implsch.o wwm_implsch_local.o wwm_implsch2.o \
       wwm_output.o wwm_compute.o wwm_diffrac.o wwm_coupl_roms.o wwm_jacobi.o wwm_coupl_shyfem.o wwm_coupl_selfe.o wwm_coupl_timor.o wwm_main.o
endif

ifdef USE_PETSC
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DPETSC
  PETSCMOD = ${PETSC_FC_INCLUDES}
  PETSCLIBS = ${PETSC_LIB}
else
  PETSCLIBS =
  PETSCMOD =
endif

ifdef USE_SED2D
  PPFLAGS := $(PPFLAGS) $(IBMCPP)-DUSE_SED2D
  OBJS := $(OBJS) sed2d_mod.o sed2d_read.o sed2d_filter.o sed2d_transport.o \
                  sed2d_friction.o sed2d_morpho.o sed2d_check.o sed2d_init.o \
                  sed2d_main.o sed2d_main_mcml.o
endif

OBJS := $(OBJS) grid_subs.o hydraulic_structures.o selfe_init.o \
        selfe_step.o selfe_finalize.o bktrk_subs.o solver_subs.o \
        io_subs.o misc_subs.o transport_TVD.o transport_TVD_imp.o sflux_9c.o lap.o

# Names of the object files with configuration intermediate file (.o,.mod) output directory prepended
CONFIG_OBJS :=  $(addprefix o/$(config)/, $(OBJS))

# Names of source files with relative path as found on $VPATH
FULLSRC := $(foreach f,$(subst .o,.F90,$(OBJS)),$(firstword $(wildcard $(VPATH:%=%$f) not-found:$f)))


################################################################################
# Library list
################################################################################
LIBS = $(MTSLIBS) $(CDFLIBS) $(GTMLIBS) $(PETSCLIBS)

################################################################################
# Build commands for SELFE
################################################################################

.PRECIOUS: $(SELFE_HOME)/mk/.check.$(config) o/$(config)/%.o o/$(config)/.dummy



# Allow for the override of executable name by the user
ifeq ($(EXEC),)
    EXEC := pelfe.$(config).ex
endif


.PHONY: pelfe

pelfe: $(EXEC) o/$(config)/.dummy



.PHONY: all
all: pelfe utility


PARMETISDIR := ./ParMetis-3.1-Sep2010

# This is the main executable target
$(EXEC): $(CONFIG_OBJS) $(PARMETISDIR)/parmetis.a $(PARMETISDIR)/metis.a
	$(FLD) $(FLDFLAGS) $(CONFIG_OBJS) $(LIBS) ${MODULEOUTDIR}. -I. -o $@


$(PARMETISDIR)/%.a:
	+$(MAKE) --directory=$(PARMETISDIR)

.PHONY: parmetis
parmetis:
	+$(MAKE) clean --directory=$(PARMETISDIR)
	+$(MAKE) --directory=$(PARMETISDIR)


.PHONY: depend
depend: o/$(config)/.dummy o/$(config)/.depend ../bin/.dummy
	@echo Dependencies generated



o/$(config)/.depend: o/$(config)/.dummy
	@echo Building depend list
	@echo $(util_srcs_abs)
	@ ../mk/sfmakedepend.pl -d o/$(config) -I $(VPATH) -f o/$(config)/.depend_extended $(FULLSRC) > depend.out
	@echo $(FULL_CONFIG_OBJS)
	@ ../mk/cull_depends.py o/$(config)/.depend_extended o/$(config)/.depend "$(FULL_CONFIG_OBJS)" 

ifneq "$(MAKECMDGOALS)" "clean"
  -include o/$(config)/.depend
endif

#o/$(config)/%.o: Hydro/%.F90 o/$(config)/.dummy o/$(config)/.depend
#	$(QUIET)$(FCP) -c $(FCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTMMOD) $(CDFMOD) $< -o $@

o/$(config)/%.o o/$(config)/%.mod: %.F90 o/$(config)/.dummy o/$(config)/.depend	
	$(QUIET)$(FCP) -c $(FCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTMMOD) $(CDFMOD) $(PETSCMOD) $< -o $(patsubst %.mod, %.o, $@)
	@touch $(patsubst %.o, %.mod, $@)




################# Utility
# The utility directory contains a lot of single-file executables plus a few that rely on shared subroutines
# for things like geometry calculations.
#
# The executables will be discovered automatically
# Files with re-usable routines should be maually added to the list util_lib_src a few lines below.
# Please consider moving the common files to UtilLib if you think the subroutines can be re-used 
# by utility scripts in several directories ... we'd like to start creating a modern API-driven library.


src_extension := .f90
util_lib_src := selfe_geometry.f90 extract_mod.f90 compute_zcor.f90 stat_pool.f90 stat_pool.f90 stripsearch_unstr.f90
util_lib_obj := $(subst $(src_extension),.o,$(util_lib_src))
util_lib_config_objs := $(addprefix o/$(config)/, $(util_lib_obj))


util_dirs :=  Combining_Scripts Grid_Scripts NCOM OneWayNestScripts Particle_Tracking Post-Processing-Fortran Sflux_nc SMS Stat3 UtilLib ArcGIS
util_dirs := $(addprefix Utility/,$(util_dirs))
VPATH := $(VPATH) $(util_dirs)


# list of source files in all the other source directories
ifneq ($(util_dirs),)
  util_srcs_abs := $(wildcard $(foreach ext,$(src_extension),$(addsuffix /*$(ext),$(util_dirs))))
  util_srcs_all := $(notdir $(util_srcs_abs) )
endif

# remove library items that don't compile as executables
util_srcs  := $(filter-out $(util_lib_src),$(util_srcs_all))
util_execs := $(addprefix ../bin/,$(subst $(src_extension),.ex,$(util_srcs)))
util_objs  := $(subst $(src_extension),.o,$(util_srcs))
util_config_objs := $(addprefix o/$(config)/, $(util_objs)) 

FULLSRC := $(FULLSRC) $(util_srcs_abs)
FULL_CONFIG_OBJS = $(CONFIG_OBJS) $(util_config_objs) $(util_lib_config_objs)

.PHONY: utility

utility: $(util_execs) o/$(config)/libutil_lib.a


#$(util_execs): $(util_config_objs) o/$(config)/libutil_lib.a
#	$(FCS) $(SCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Lo/$(config) -lutil_lib -Io/$(config)  $(CDFLIBS) $(GTTMOD) $(CDFMOD) $< -o $@
$(util_execs): ../bin/%.ex : o/$(config)/%.o o/$(config)/libutil_lib.a
	$(FCS) $(SCPFLAGS) $< ${MODULEOUTDIR}o/$(config) -Lo/$(config) -lutil_lib -Io/$(config)  $(CDFLIBS) $(GTTMOD) $(CDFMOD) -o $@

o/$(config)/libutil_lib.a: $(util_lib_config_objs)
	$(AR) $(ARFLAGS) $@ $^


#$(util_config_obs) : $(util_srcs) o/$(config)/.dummy  ../bin/.dummy
#	$(SCP) -c $(FCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTTMOD) $(CDFMOD) $< -o $@


o/$(config)/%.o : %.f90 o/$(config)/.dummy  ../bin/.dummy
#	$(FCS) -c $(FCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTTMOD) $(CDFMOD) $< -o $@
	$(FCS) -c $(FCPFLAGS) $< ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTTMOD) $(CDFMOD) -o $@


../bin/%.ex: o/$(config)/%.o  ../bin/.dummy
	$(FLD) $(FLDFLAGS) $< $(LIBS) ${MODULEOUTDIR}. -I. -o $@





#################################################################################

# This file exists solely to force creation of the subdirectories that
# hold the dependency, object and preprocessor files.
# NOTE: we cant use a dependency for this because .dummy will always
#       be out-of-date with respect to the subdirs because they are
#       modified after .dummy is created.  But since there are no
#       dependencies, this rule only runs when the file doesn't exist.
outdir_exists := $(wildcard o/$config)
make_defs_local_exists := $(wildcard ../mk/Make.defs.local)

o/$(config)/.dummy:
	@echo $(make_defs_local_exists)
        ifeq ( $(strip $(make_defs_local_exists)), )
             $(error Make.defs.local does not exist in the /mk directory. Please create one using another mk/Make.defs.xxx as a template.  \
             If you put your Make.defs.xxx file in svn, you can make a symbolic link Make.defs.local pointing to it.[1])
        endif
	$(ECHO)mkdir -p o/$(config)
	@touch o/$(config)/.dummy


../bin/.dummy:
	$(ECHO)mkdir -p ../bin
	@touch ../bin/.dummy

../mk/Make.defs.local :
	$(error Make.defs.local does not exist in the /mk directory. Please create one using another mk/Make.defs.xxx as a template. \n \
             If you put your Make.defs.xxx file in svn, you can make a symbolic link Make.defs.local pointing to it.[1])



################################################################################
# Clean up or Clobber
################################################################################
# Clean up objects and modules. Clean kills all configurations.

.PHONY: clean
clean:
	$(QUIET)rm -f *.F90 *.c *.o *.mod $(EXEC)
	$(QUIET)rm -rf o/*
	$(QUIET)rm -rf ../bin/*
	$(QUIET)rm -rf $(PARMETISDIR)/*.a

