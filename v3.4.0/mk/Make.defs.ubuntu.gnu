################################################################################
# Parallel SELFE Makefile
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


ENV = ubuntu.gnu

################################################################################
# Environment for Ubuntu with GNU compilers
# Packages requirement: gfortran, gcc, libnetcdf-dev, mpich2, libparmetis-dev, libcr-dev,
# which can be installed with the command below:
# "sudo apt-get install gfortran gcc libnetcdf-dev mpich2 libparmetis-dev libcr-dev"
# Test is passed on Ubuntu 11.10 64bit, Ubuntu 12.04 64bit
################################################################################

#  USE_WRAP = yes
FCP = /usr/bin/mpif90.mpich2 -f90=gfortran -ffree-line-length-none
FLD = $(FCP)
# MPI vserion (1 or 2) 
PPFLAGS := $(PPFLAGS) -DMPIVERSION=2 #-DUSE_WRAP

FCPFLAGS = $(PPFLAGS) -O2 -Bstatic #-g -fbacktrace #-finit-real=nan -fbounds-check #MPI code
FLDFLAGS = -O2 #for final linking of object files

  #####Libraries
#MTSLIBS = -L/usr/lib/ -lparmetis -lmetis
MTSLIBS = -L./ParMetis-3.1-Sep2010 -lparmetis -lmetis
CDFLIBS = -L/usr/lib/ -lnetcdf -lnetcdff
CDFMOD = -I/usr/include/ # modules for netcdf
ifdef USE_GOTM
   GTMMOD =  -I/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/modules/IFORT/ #modules
   GTMLIBS = -L/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/lib/IFORT/ -lturbulence_prod  -lutil_prod
else
   GTMMOD =
   GTMLIBS =
endif



################################################################################
# Alternate executable name if you do not want the default. 
################################################################################

#EXEC   := othername.ex


################################################################################
# Algorithm preference flags.
# Comment out unwanted modules and flags.
################################################################################

# -DSELFE is always on and is defined elsewhere

# Precip/evaporation model
# PPFLAGS := $(PPFLAGS) -DPREC_EVAP 

# MM5 in heat exchange model
# PPFLAGS := $(PPFLAGS) -DMM5

# TVD flag (turn off for performance)
#  CHOOSE_TVD = yes
#   EXEC := $(EXEC)_CHOOSE_TVD

# GOTM turbulence closure model
# USE_GOTM = yes

# Wind wave model WWM
# USE_WWM = yes

# Implicit schemes in WWM needs PetSc
# USE_PETSC = yes

# TIMOR 
# USE_TIMOR = yes

# Harmonic analysis tool
# USE_HA = yes

##### Select only _one_ model from below

# Ecological model - NAPZD Spitz (internal use only)
# USE_NAPZD = yes

# Or:
# Ecological model (EcoSim)
# USE_ECO = yes

# Or:
# CE-QUAL-ICM
# USE_ICM = yes

# Or:
# Sediment model 
# USE_SED = yes

# If you choose SED you should set the following algorithmic preferences

  ##Bedload 
#  PPFLAGS := $(PPFLAGS) -DBEDLOAD

  ##Bedload - MPM model
#  PPFLAGS := $(PPFLAGS) -DBEDLOAD_MPM

  ##slope formulation
#  PPFLAGS := $(PPFLAGS) -DDAMGAARD
#  PPFLAGS := $(PPFLAGS) -DDELFT
#  PPFLAGS := $(PPFLAGS) -DCARMO

  ##Bedload - VR model
#  PPFLAGS := $(PPFLAGS) -DBEDLOAD_VR

  ##Suspended load
#  PPFLAGS := $(PPFLAGS) -DSUSPLOAD

  ##boundary conditions for WENO
  ## default strictly monotonic
#  PPFLAGS:= $(PPFLAGS) -DLINEAR_CONTINUATION
#  PPFLAGS:= $(PPFLAGS) -DNEUMANN

  ## Morphology
#  PPFLAGS := $(PPFLAGS) -DSED_MORPH

  ## Choose one drag formulation from the following 3 choices (only 1st one is functional now)
#  PPFLAGS := $(PPFLAGS) -DUV_LOGDRAG
  #PPFLAGS := $(PPFLAGS) -DUV_QDRAG
  #PPFLAGS := $(PPFLAGS) -DUV_LDRAG

  ##sediment density in eqstate
#   PPFLAGS:= $(PPFLAGS) -DDENSED

# Or:
# Sediment model (2D)
# #  USE_SED2D = yes

# Or:
# Oil spill model (not active)
# USE_OIL = yes


#########  Compiler configuration related flags

# Include a timer
# TIMER = yes

# Debug mode (more time consuming)
# USE_DEBUG = yes


######### Specialty compiler flags and workarounds
# Add -DNO_TR_15581 like below for allocatable array problem in sflux_subs.F90
# PPFLAGS := $(PPFLAGS) -DNO_TR_15581

# Obsolete flags: use USE_WRAP flag to avoid problems in ParMetis lib (calling C from FORTRAN)
# PPFLAGS := $(PPFLAGS) -DUSE_WRAP 



#############################################################################################

