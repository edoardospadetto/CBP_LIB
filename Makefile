CC = gcc
CFOR = gfortran
CFORMPI = mpif90
SOURCES=gemm.cpp
DFLAGS = -ggdb3
LIBSPATH=-L/cineca/prod/opt/compilers/cuda/10.1/none/update2/lib64
LIBS = $(LIBSPATH) -lcuda -lcudart -lcublas
EXEC = -o gemm
TYPE = -D CBP$(NAME)
all :
	#$(CC) $(SOURCES) $(LIBS)  $(EXEC)


	$(CC) $(DFLAGS) -c  ./src/cbp_wrap.cpp -o ./build/cbp_wrap.o $(TYPE) -lblas
	$(CFOR) $(DFLAGS) -c ./src/cbp_parallel.F90 -o ./build/cbp_parallel.o $(TYPE)
	$(CFORMPI) $(DFLAGS) -c ./src/cbp_distribute_mod.F90 -o ./build/cbp_distribute_mod.o $(TYPE)
	#$(CFORMPI) $(DFLAGS) -c ./test_mpi_perform.F90 -o ./build/test_mpi_perform.o $(TYPE)
	$(CFORMPI) $(DFLAGS) -c ./basic_test.F90 -o ./build/basic_test.o $(TYPE)	#Works only with real and double
	$(CFORMPI)  $(DFLAGS) $(LIBGFORTR) ./build/basic_test.o ./build/cbp_parallel.o ./build/cbp_wrap.o ./build/cbp_distribute_mod.o -lstdc++  $(LIBS)   -o test$(NAME)
	#$(CFORMPI)  $(DFLAGS) $(LIBGFORTR) ./build/test_mpi_perform.o ./build/cbp_parallel.o ./build/cbp_wrap.o ./build/cbp_distribute_mod.o -lstdc++  $(LIBS)   -o ./bin/test$(NAME)
	rm ./build/*

