#!/usr/bin/bash

ulimit -s unlimited

#EXPORTING OPENMP ENVIRONMENT VARIABLES
export OMP_NUM_THREADS=20
export MKL_NUM_THREADS=20
export OMP_SCHEDULE="dynamic"
# export OMP_STACKSIZE=16M
# export OMP_PROC_BIND="close"


# -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all \
# -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan \
# -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all \
# -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan \-fno-trapping-math -fno-signed-zeros 

#COMPILING SOURCE FILES#

#GNU COMPILERS

gfortran-12 -ffree-line-length-512 -Ofast -fopenmp -march=native -mtune=native \
-mavx2 -c -funroll-loops -flto=auto -ftree-vectorize -fstack-arrays -finline-functions \
-m64  -I"${MKLROOT}/include" \
initialize.f90

gfortran-12 -ffree-line-length-512 -Ofast -fopenmp -march=native -mtune=native \
-mavx2 -std=f2008 -c -funroll-loops -flto=auto -ftree-vectorize -fstack-arrays -finline-functions \
particle.f90 interactions.f90 transfer.f90 domain.f90 kernel.f90 functions.f90 \
search.f90 gradient.f90 setup.f90 solver.f90 boundary.f90 output.f90 internal.f90 \
memory.f90 isph.f90 integrator.f90 turbulence.f90 part_shift.f90 \
scalar.f90 porous.f90

# ifx -Ofast -qopenmp -xHost -std08 -ipo -c -qmkl=parallel -heap-arrays \
# particle.f90 interactions.f90 transfer.f90 domain.f90 kernel.f90 functions.f90 \
# search.f90 gradient.f90 setup.f90 solver.f90 boundary.f90 output.f90 internal.f90 \
# memory.f90 isph.f90 integrator.f90 turbulence.f90 part_shift.f90 \
# mls.f90 scalar.f90 porous.f90
# ifx  -Ofast -qopenmp -xHost -ipo -c -qmkl=parallel -heap-arrays \
# initialize.f90

#COMPILING EXECUTABLE#

#GNU COMPILERS

# gfortran-12 -ffree-line-length-512 -Ofast -fprefetch-loop-arrays -fopenmp -march=native -mtune=native \
# -mavx2 -std=f2008 -funroll-loops -flto=auto -ftree-vectorize -fstack-arrays \
# -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_intel_thread \
# -lmkl_core -liomp5 -lpthread -lm -ldl -fwhole-program \
# particle.o interactions.o transfer.o domain.o initialize.o kernel.o functions.o \
# search.o gradient.o setup.o solver.o boundary.o output.o internal.o \
# memory.o isph.o integrator.o turbulence.o part_shift.o \
# mls.o scalar.o porous.o main.f90 -o main


gfortran-12 -ffree-line-length-512 -Ofast -fopenmp -march=native -mtune=native \
-mavx2 -std=f2008 -funroll-loops -flto=auto -ftree-vectorize -fstack-arrays -finline-functions \
particle.o interactions.o transfer.o domain.o initialize.o kernel.o functions.o \
search.o gradient.o setup.o solver.o boundary.o output.o internal.o \
memory.o isph.o integrator.o turbulence.o part_shift.o \
scalar.o porous.o main.f90 -o main \
-Wl,--start-group ${MKLROOT}/lib/libmkl_gf_lp64.a ${MKLROOT}/lib/libmkl_intel_thread.a \
${MKLROOT}/lib/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl


# ifx  -Ofast -qopenmp -xHost -std08 -ipo -qmkl=parallel -heap-arrays \
# particle.o interactions.o transfer.o domain.o initialize.o kernel.o functions.o \
# search.o gradient.o setup.o solver.o boundary.o output.o internal.o \
# memory.o isph.o integrator.o turbulence.o part_shift.o \
# mls.o scalar.o porous.o main.f90 -o main -lm -ldl
