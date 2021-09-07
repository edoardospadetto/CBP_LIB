CFOR = mpif90

FILES = misc_mod debug_module matrix_interface scalapack_interface $(TEST)_test


SRC=./src/
BUILD= ./build/
LIBS =  -L $(SCALAPACK_LIB) -L $(LAPACK_LIB) -L$(BLAS_LIB) -lscalapack -llapack -lblas
DEF = -D CBP$(TYPE) -O3


all :  $(patsubst %, $(BUILD)%.o, $(FILES)) 
	$(CFOR) $(DEF)  $^ $(LIBS)  -o $(TEST)_exec

$(BUILD)%.o: $(SRC)%.F90
	$(CFOR)  $(DEF) $(LIBS) -c $(CFLAGS) $< -o $@
	


clean: 
	rm *.mod
	rm $(BUILD)*.o
	rm *_exec
