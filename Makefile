# Makefile to build interpinic on various platforms
# Note: If netcdf library is not built in the standard location, 
# you must set the environment variables INC_NETCDF and LIB_NETCDF

EXEDIR = .
EXENAME = richards
RM = rm

.SUFFIXES:
.SUFFIXES: .f90 .o

# Determine platform 
UNAMES := $(shell uname -s)
UNAMEM := $(findstring CRAY,$(shell uname -m))

# Architecture-specific flags and rules
#
#------------------------------------------------------------------------
# Cray 
#------------------------------------------------------------------------

ifeq ($(UNAMEM),CRAY)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SGI
#------------------------------------------------------------------------

ifeq ($(UNAMES),IRIX64)
FC = f90
FFLAGS = -64 -c -trapuv -I$(INC_NETCDF) -g -C -DEBUG:trap_uninitialized=ON
LDFLAGS = -64 -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# SUN
#------------------------------------------------------------------------

ifeq ($(UNAMES),SunOS)
FC = f90
FFLAGS = -c -stackvar -f -I$(INC_NETCDF) -g
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# AIX
#------------------------------------------------------------------------

ifeq ($(UNAMES),AIX)
LIB_NETCDF := /usr/local/apps/netcdf-3.5/lib32/r4i4
FC = xlf90
FFLAGS = -c -I$(INC_NETCDF) -qsuffix=f=f90 -g -qfullpath
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# OSF1
#------------------------------------------------------------------------

ifeq ($(UNAMES),OSF1)
FC = f90
FFLAGS = -c -I$(INC_NETCDF)
LDFLAGS = -L$(LIB_NETCDF) -lnetcdf
endif

#------------------------------------------------------------------------
# Linux
#------------------------------------------------------------------------

ifeq ($(UNAMES),Linux)
FC = gfortran
FFLAGS =  -c -g
LDFLAGS = -g
endif

#------------------------------------------------------------------------
# Default rules and macros
#------------------------------------------------------------------------

OBJS := richardsmain.o saxtonscheme.o soitexmod.o tridiagonalmod.o iomod.o

.f90.o:
	$(FC) $(FFLAGS) $<

$(EXEDIR)/$(EXENAME): $(OBJS)
	$(FC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	$(RM) -f $(OBJS) *.mod $(EXEDIR)/$(EXENAME)

richardsmain.o : richardsmain.f90 soitexmod.o saxtonscheme.o tridiagonalmod.o iomod.o
saxtonscheme.o : saxtonscheme.f90 soitexmod.o
soitexmod.o : soitexmod.f90
tridiagonalmod.o : tridiagonalmod.f90
iomod.o : iomod.f90
