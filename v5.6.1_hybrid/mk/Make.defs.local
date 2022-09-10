################################################################################
# Parallel SCHISM Makefile
#
# User makes environment settings for particular OS / PLATFORM / COMPILER / MPI
# below as well as setting flags having to do with included algorithms (e.g. sediment)
# and the compiler configuration (debug, timing). 
#
# The environment settings are based on the following options.
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
# Libraries (needed for parallel code)
#   MTSLIBS: Flags for linking ParMeTiS/MeTiS libaries
#   ALTLIBS: Flags for linking alternate solver libraries (LAPACK or ITPACK,
#            these are just for testing)
#
#
#
################################################################################
ENV = thor.intel

################################################################################
# Environment for Linux / 64 bit /  Intel Compiler / MPICH2 (Ymir cluster)
################################################################################

#FCP = mpif90 -f90=ifort
FCP = mpiifort -f90=ifort
FLD = $(FCP)
# MPI version (1 or 2)
PPFLAGS := $(PPFLAGS) -DMPIVERSION=2 #-DUSE_WRAP
OPTFLAGS = -O2
FCPFLAGS = $(PPFLAGS) $(OPTFLAGS) #-g -fbacktrace #-finit-real=nan -fbounds-check #MPI code
FLDFLAGS = $(OPTFLAGS) #for final linking of object files
#PPFLAGS := $(PPFLAGS) -DMPIVERSION=2
#FCPFLAGS = $(PPFLAGS) -O1 -g -fp-model precise -traceback -assume byterecl -check uninit -check bounds -check pointers -warn all,nounused -gen-interfaces 
#FLDFLAGS = -O1 -g -fp-model precise -traceback -assume byterecl -check uninit -check bounds -check pointers -warn interfaces,nouncalled -gen-interfaces 

# Libraries
# Parmetis
MTSLIBS = -L./ParMetis-3.1-Sep2010 -lparmetis -lmetis
# Netcdf installed locally
# export NCDIR=/home/kmartins/DEV/Libs/netcdfC
# export NFDIR=/home/kmartins/DEV/Libs/netcdfF
# Netcdf installed on Thor
export NCDIR=/sw/lib/netcdf/ifort/4.4.4_linked_with_4.5.0
export NFDIR=/sw/lib/netcdf/ifort/4.4.4_linked_with_4.5.0
CDFLIBS=-L${NCDIR}/lib -lnetcdff -lnetcdf
CDFMOD=-I${NFDIR}/include

#ifdef USE_GOTM
#   GTMMOD =  -I/home/tgueri02/SCHISM_v5_3/v5.3/src/GOTM3.2.5/netcdf_include #modules for netcdf
#   GTMLIBS = -L/Utilisateurs/kli/opt/source/gotm-4.0.0/lib/IFORT -lturbulence_prod  -lutil_prod
#else
#   GTMMOD =
#   GTMLIBS =
#endif


################################################################################
# Alternate executable name if you do not want the default. 
################################################################################

EXEC   := selfe


################################################################################
# Algorithm preference flags.
# Comment out unwanted modules and flags.
################################################################################

# -DSCHISM is always on and is defined elsewhere
include ../mk/include_modules

# Don't comment out the follow ifdef
# Note: currently GOTM4 may give reasonable results only with k-omega
ifdef USE_GOTM
  #Following for GOTM4
  #GTMMOD =  -I/sciclone/home04/yinglong/SELFE/svn/trunk/src/GOTM4.0/modules/PGF90/ #modules
  #GTMLIBS = -L/sciclone/home04/yinglong/SELFE/svn/trunk/src/GOTM4.0/lib/PGF90/ -lturbulence_prod -lutil_prod

  #Following for GOTM3
  GTMMOD =  -I/sciclone/home04/yinglong/gotm-3.2.5/modules/PGF90/ #modules
  GTMLIBS = -L/sciclone/home04/yinglong/gotm-3.2.5/lib/PGF90/ -lturbulence_prod -lutil_prod
else
  GTMMOD =
  GTMLIBS =
endif


######### Specialty compiler flags and workarounds
# Add -DNO_TR_15581 like below for allocatable array problem in sflux_subs.F90
# PPFLAGS := $(PPFLAGS) -DNO_TR_15581

# Obsolete flags: use USE_WRAP flag to avoid problems in ParMetis lib (calling C from FORTRAN)
# PPFLAGS := $(PPFLAGS) -DUSE_WRAP 



#############################################################################################

