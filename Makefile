CC=cilk++
CFLAGS=-O3 

all: boruvka hybrid hybrid_omp

boruvka:
	@echo "Generating executables for cuda implementation"
	nvcc -I$(HOME)/cudpp/include -arch=compute_35 -code=sm_35 boruvka.cu -o boruvka -lcudpp -L$(HOME)/cudpp/build/lib
hybrid:
	@echo "Generating executables for our hybrid approach"
	nvcc -I$(HOME)/cudpp/include -arch=compute_35 -code=sm_35 hybrid.cu -o hybrid -lcudpp -L$(HOME)/cudpp/build/lib
hybrid_omp:
	@echo "Generating executables for hybrid with omp implementation"
	nvcc -I$(HOME)/cudpp/include -arch=compute_35 -code=sm_35 -Xcompiler -fopenmp hybrid_with_omp.cu -o hybrid_omp -lgomp -lcudpp -L$(HOME)/cudpp/build/lib
clean:
	@echo "Cleaning up all the generated executables"
	rm hybrid boruvka hybrid_omp
