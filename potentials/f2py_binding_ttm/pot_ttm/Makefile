FC=gfortran
FPY=f2py

SRC=$(addprefix ./, \
	model_mod.f90 \
	ttm2f_mod.f90 \
	ttm3f_mod.f90 \
	nasa_mod.f90 \
	math.f90 \
	smear.f90 \
	read_conf.f90 \
	nasa.f90 \
	potential_mod.f90 \
	dfpmin.f90)

FPYSRC=$(addprefix ./, \
	ttm_from_f2py.f90)

OBJ=$(addprefix ./, \
	model_mod.o \
	ttm2f_mod.o \
	ttm3f_mod.o \
	nasa_mod.o \
	math.o \
	smear.o \
	read_conf.o \
	nasa.o \
	potential_mod.o \
	dfpmin.o )

SO=./ttm_all.so
MODULE=ttm

FFLAGS=-O2 -fPIC
LINKFLAGS= -O2 -fPIC -shared
FPYCOMMAND= -c $(FPYSRC) -m $(MODULE) $(SO)

library: $(OBJ)
		$(FC) -o $(SO) $(LINKFLAGS) $(OBJ)
all: $(OBJ)
			$(FPY) $(FPYCOMMAND)
$(OBJ) : Makefile

%.o:%.f90
	        $(FC) $(FFLAGS) -c $< -o $@
%.o:%.F90
	        $(FC) $(FFLAGS) -c $< -o $@

clean: 
	rm -rf $(OBJ) *.so *.o *.mod
