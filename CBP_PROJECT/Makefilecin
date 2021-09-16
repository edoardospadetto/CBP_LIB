CC = gcc
CFOR =  pgf90 -x=f90  -Mpreprocess -Mcuda -gpu=rdc
CFORMPI = mpif90 -cpp  -Mcuda
SOURCES=gemm.cpp
DFLAGS =  -g
LIBSPATH= -I /home/edoardo/lib/Linux_x86_64/21.7/cuda/include -I /home/edoardo/lib/Linux_x86_64/21.7/math_libs/10.2/include -L /home/edoardo/lib/Linux_x86_64/21.7/math_libs/10.2/lib64
#-L/cineca/prod/opt/compilers/cuda/10.1/none/update2/lib64
LIBS = $(LIBSPATH) -lcublas -lcuda -lcudart  -lcusparse -llapack
EXEC = -o gemm
TYPE = -D CBP$(NAME)
SRC = ./src/CBP_lib

all :  cbp_parallel distribute

unittest : clean distribute
	$(CFORMPI) $(DFLAGS) $(LIBSPATH) -c ./main_s/test_mpi_perform.F90 -o ./build/test_mpi_perform.o $(TYPE)
	$(CFORMPI)  $(DFLAGS) $(LIBGFORTR) ./build/test_mpi_perform.o ./build/cbp_parallel.o ./build/cbp_wrap.o ./build/cbp_distribute_mod.o -lstdc++  $(LIBS)   -o ./bin/test$(NAME)

lapack :
			$(CFOR)  $(TYPE) $(DFLAGS) $(LIBS) -c $(CFLAGS) src/lapack/src/lapack_interface.f90 -o ./build/lapack_interface.o

kernels :

	echo 'Compiling Kernel..'

	pgf90 -x=f90  -Mpreprocess -Mcuda -gpu=rdc -c $(TYPE) $(SRC)/kernels.cuf -o ./build/kernels.o

easytest : all misc lapack
	echo ' '
	echo Linking..
	echo ' '
	$(CFORMPI) $(DFLAGS) -c ./main_s/$(WHICH)_test.F90 -o ./build/$(WHICH)_test.o $(TYPE)
	$(CFORMPI)  $(DFLAGS) $(LIBGFORTR)  ./build/$(WHICH)_test.o ./build/misc_mod.o ./build/cbp_parallel.o ./build/cbp_wrap.o ./build/lapack_interface.o ./build/cbp_distribute_mod.o ./build/kernels.o -lstdc++  $(LIBS)   -o ./bin/test$(WHICH)$(NAME)$(PROCS)


temporal : all misc
	echo ' '
	echo Linking..
	echo ' '
	$(CFORMPI) $(DFLAGS) -c ./temporal_test/$(WHICH)_time.F90 -o ./build/$(WHICH)_time.o $(TYPE)
	$(CFORMPI)  $(DFLAGS) $(LIBGFORTR)   ./build/$(WHICH)_time.o ./build/misc_mod.o ./build/cbp_parallel.o ./build/cbp_wrap.o ./build/cbp_distribute_mod.o ./build/kernels.o -lstdc++  $(LIBS)   -o ./bin/testb$(WHICH)$(NAME)$(PROCS)



distribute : kernels cbp_parallel
	$(CFORMPI) $(DFLAGS) -c $(SRC)/cbp_distribute_mod.F90 -o ./build/cbp_distribute_mod.o $(TYPE)

cbp_parallel:
	$(CFOR) $(DFLAGS) -c $(SRC)/cbp_parallel.F90 $(LIBS) -o ./build/cbp_parallel.o $(TYPE)
	$(CC) $(DFLAGS)  $(LIBSPATH) -c $(SRC)/cbp_wrap.cpp -o ./build/cbp_wrap.o $(TYPE) -lblas

misc:
	$(CFOR) $(DFLAGS) -c ./src/misc_mod.f90 $(LIBS) -o ./build/misc_mod.o $(TYPE)



clean : cleanbin
	rm ./build/*
	rm *.mod


cleanbin:
	rm ./bin/*
