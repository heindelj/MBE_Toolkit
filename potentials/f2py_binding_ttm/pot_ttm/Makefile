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

#EXE=main.x
SO=./ttm_all.so
MODULE=ttm

FFLAGS=-O2 -shared -fPIC
PYFLAGS=-c -m

library: $(OBJ)
	        $(FC) $(FFLAGS) $(OBJ) $(INCLUDE_PATH) $(LIB_PATH) $(LIBSHDF) $(LIBS) -o $(SO)

all: $(OBJ)
			$(FPY) $(PYFLAGS) $(MODULE) $(FPYSRC) $(SO)

$(OBJ) : Makefile

%.o:%.f90
	        $(FC) $(FFLAGS) $(DEFS) $(INCLUDE_PATH) -c $< -o $@
%.o:%.F90
	        $(FC) $(FFLAGS) $(DEFS) $(INCLUDE_PATH) -c $< -o $@

clean: 
	rm -rf $(OBJ) *.so *.o *.mod
