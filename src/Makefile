CC = gcc
CFOR = gfortran
SOURCES=gemm.cpp
DFLAGS = -ggdb3
LIBSPATH=-L/cineca/prod/opt/compilers/cuda/10.1/none/update2/lib64
LIBS = $(LIBSPATH) -lcuda -lcudart -lcublas
TYPE = -D CBPREAL
EXEC = -o gemm
all :
	#$(CC) $(SOURCES) $(LIBS)  $(EXEC)


	$(CC) $(DFLAGS) -c  ./cbp_wrap.cpp -o ./cbp_wrap.o $(TYPE)
	$(CFOR) $(DFLAGS) -cpp -c ./gemm_example.F90 -o ./gemm_example.o $(TYPE)
	$(CFOR)  $(DFLAGS) $(LIBGFORTR) ./gemm_example.o ./cbp_wrap.o -lstdc++  $(LIBS)   -o test

