############################################################################
#  QUEST Makefile    
#    Date:   08/03/2014
#  Author:   Chia-Chen Chang
#
#  Modified for Apple Silicon (arm64) macOS with Homebrew:
#    - removed -m64 (not accepted by gcc on arm64)
#    - added -fallow-argument-mismatch (needed by gfortran >= 10 for
#      legacy BLAS/LAPACK calls in QUEST)
#    - C++ compiler taken from Homebrew gcc (g++-NN), so that -lstdc++
#      resolves when linking with gfortran
#    - new LAPACK option "brew": Homebrew OpenBLAS instead of the bundled
#      (old, x86-only) OpenBLAS
############################################################################
QUEST_DIR = $(shell pwd)

# Homebrew prefix (/opt/homebrew on Apple Silicon)
BREW_PREFIX = $(shell brew --prefix 2>/dev/null)

# 1) gnu, 2) intel
COMPILER  = gnu

# 1) default, 2) mkl_seq, 3) mkl_par, 4) brew (Homebrew OpenBLAS, use this on a Mac)
LAPACK    = brew

# intel MKL library path
MKLPATH   = $(MKLROOT)/lib/intel64

# MAGMA library path
MAGMAPATH = 

# nVidia CUDA installation path
CUDAPATH  = 

# Checkboard decomposition
FLAG_CKB  = #-DDQMC_CKB

# GPU version equal-time Green's function kernel
FLAG_ASQRD = -DDQMC_ASQRD

# GPU version time-dependent Green's function kernel
FLAG_BSOFI = #-DDQMC_BSOFI

# Enabling nVidia CUDA support in DQMC
FLAG_CUDA = #-DDQMC_CUDA



# --------------------------------------------------------------------------
#  Check ASQRD and CKB compatibility
# --------------------------------------------------------------------------
ifdef FLAG_ASQRD
  ifdef FLAG_CKB
    $(error ASQRD method does not support checkerboard decomposition at the moment. \
      Please turn off the flag FLAG_CKB or FLAG_ASQRD)
  endif
endif



# --------------------------------------------------------------------------
#  Compiler and flags
# --------------------------------------------------------------------------
ifndef COMPILER
  $(error "COMPILER" is not defined in the Makefile.)
endif

ifeq ($(COMPILER), intel)
  FC        = ifort
  CXX       = icpc
  FC_FLAGS  = -m64 -warn all -O3 -unroll
  #FC_FLAGS = -m64 -g -traceback -check all -O0 -ftrapuv -debug all
  #CXX_FLAGS = -m64 -g -traceback -O0 -check-uninit -ftrapuv -debug all
  CXX_FLAGS = -m64 -Wall -O3 -unroll $(CUDAINC) $(MAGMAINC)
endif
ifeq ($(COMPILER), gnu)
  FC        = gfortran
  # Homebrew's GNU g++ (e.g. g++-14); plain "g++" on a Mac is Apple clang
  CXX       = $(lastword $(sort $(wildcard $(BREW_PREFIX)/bin/g++-[0-9]*)))
  ifeq ($(CXX),)
    $(error Homebrew g++ not found. Run "brew install gcc")
  endif
  FC_FLAGS  = -Wall -Wno-tabs -O3 -funroll-loops -fallow-argument-mismatch
  CXX_FLAGS = -Wall -O3 -funroll-loops $(CUDAINC) $(MAGMAINC)
endif



# --------------------------------------------------------------------------
# C++ libraries
# --------------------------------------------------------------------------
CXXLIB = -lstdc++ #-lrt



# --------------------------------------------------------------------------
#  BLAS and LAPACK
# --------------------------------------------------------------------------
ifndef LAPACK
  $(error "LAPACK" is not defined in the Makefile.)
endif

ifeq ($(LAPACK), default)
  libOpenBLAS   = $(QUEST_DIR)/OpenBLAS/libopenblas.a

  # Pass if the target is 'clean' or 'lapack'
  ifneq ($(MAKECMDGOALS), clean)
  ifneq ($(MAKECMDGOALS), lapack)
    ifeq ($(wildcard $(libOpenBLAS)),)
      $(error It appears that LAPACK is set to be "default", but libopenblas.a is missing. \
        Use "make lapack" to build the required library.)
    endif
  endif
  endif
  LAPACKLIB = $(libOpenBLAS)
endif

ifeq ($(LAPACK), brew)
  OPENBLAS_DIR = $(shell brew --prefix openblas 2>/dev/null)
  ifneq ($(MAKECMDGOALS), clean)
    ifeq ($(wildcard $(OPENBLAS_DIR)/lib/libopenblas*),)
      $(error Homebrew OpenBLAS not found. Run "brew install openblas")
    endif
  endif
  LAPACKLIB = -L$(OPENBLAS_DIR)/lib -lopenblas
endif

ifeq ($(LAPACK), mkl_seq)
  ifdef MKLPATH
    LAPACKLIB = -L$(MKLPATH) -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
  else
    $(error "MKLPATH" is not defined in the Makefile.)
  endif
endif

ifeq ($(LAPACK), mkl_par)
  ifdef MKLPATH
    LAPACKLIB = -L$(MKLPATH) -mkl:parallel -lpthread
  else
    $(error "MKLPATH" is not defined in the Makefile.)
  endif
endif


# --------------------------------------------------------------------------
# MAGMA library path
# --------------------------------------------------------------------------
ifdef MAGMAPATH
  MAGMALIB = -L$(MAGMAPATH)/lib -lmagma -lmagmablas -lmagma 
  MAGMAINC = -I$(MAGMAPATH)/include
endif



# --------------------------------------------------------------------------
# CUDA compiler and libraries
# --------------------------------------------------------------------------
ifdef CUDAPATH
  ifeq ($(MAGMAINC),)
    $(error In Makefile, variable MAGMATH is not defined in Makefile)
  endif
  NVCC = $(CUDAPATH)/bin/nvcc
  CU_FLAGS = -O3 -Xptxas -v -m 64 -arch sm_20 $(MAGMAINC)
  CUDALIB = -L$(CUDAPATH)/lib64 -lcublas -lcudart -lcuda
  CUDAINC = -I$(CUDAPATH)/include
endif


# --------------------------------------------------------------------------
# DQMC library
# --------------------------------------------------------------------------
DQMCLIB    = libdqmc.a


# --------------------------------------------------------------------------
# Libraries required by driver routines 
# --------------------------------------------------------------------------
LIB        = $(CXXLIB) $(LAPACKLIB) $(CUDALIB) $(MAGMALIB)

INC        = $(MPIINC)


# --------------------------------------------------------------------------
# Archiver and flags
# --------------------------------------------------------------------------
ifeq ($(COMPILER), intel) 
	ARCH       = xiar
else ifeq ($(COMPILER), gnu)
	ARCH = ar
endif
ARFLAG     = cr
RANLIB     = ranlib



# --------------------------------------------------------------------------
# Optional complication flags
# -D_SXX, -D_QMC_MPI
# --------------------------------------------------------------------------
PRG_FLAGS = $(FLAG_BSOFI) $(FLAG_ASQRD) $(FLAG_CKB) $(FLAG_CUDA) 

FLAGS = $(FC_FLAGS) $(PRG_FLAGS)



# --------------------------------------------------------------------------
export 

.PHONY: libdqmc example libopenblas clean

all : libdqmc example

libdqmc :
	(cd SRC; $(MAKE))
example : libdqmc
	(cd EXAMPLE; $(MAKE))

lapack : libopenblas

libopenblas : 
	(cd $(QUEST_DIR)/OpenBLAS; $(MAKE))

clean :
	-(cd $(QUEST_DIR)/OpenBLAS; $(MAKE) clean)
	(cd $(QUEST_DIR)/SRC; $(MAKE) clean)
	(cd $(QUEST_DIR)/EXAMPLE; $(MAKE) clean)
	(rm -f $(QUEST_DIR)/$(DQMCLIB))
