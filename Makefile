STD_VERSION = -std=c++17
INC_FLAGS = -Iinc
NVCC_INC_FLAGS = -gencode=arch=compute_70,code=\"sm_70\"
LD_FLAGS = -lcudart -L/usr/local/cuda/lib64
CPP_HEADER_FILES = inc/misc_def.h inc/host_data_structures.hpp inc/device_data_structures.h inc/file_parsing.hpp inc/simulation.h
CPP_SOURCE_FILES = src/main.cpp src/file_parsing.cpp
EXE = cuAleae

test: ${CPP_HEADER_FILES} 
	nvcc -G -g ${NVCC_INC_FLAGS} ${INC_FLAGS} -c src/simulation_kernel.cu -o simulation.o
	g++ -g ${STD_VERSION} ${INC_FLAGS} ${CPP_SOURCE_FILES} simulation.o -o ${EXE} ${LD_FLAGS}

build: ${CPP_HEADER_FILES}
	nvcc -O2 ${NVCC_INC_FLAGS} ${INC_FLAGS} -c src/simulation_kernel.cu -o simulation.o
	g++ -O2 ${STD_VERSION} ${INC_FLAGS} ${CPP_SOURCE_FILES} simulation.o -o ${EXE} ${LD_FLAGS}

clean:
	rm -rf $(EXE)
	rm -rf *.o