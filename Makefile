STD_VERSION = -std=c++17
INC_FLAGS = -Iinc

# Edit these flags and add sm version compatible with your nVidia GPU
NVCC_INC_FLAGS = -arch=sm_75 \
		-gencode=arch=compute_75,code=sm_75 \
		-gencode=arch=compute_80,code=sm_80 \
		-gencode=arch=compute_86,code=sm_86 \
		-gencode=arch=compute_87,code=sm_87 \
		-gencode=arch=compute_90,code=sm_90 \
		-gencode=arch=compute_100,code=sm_100 \
		-gencode=arch=compute_120,code=sm_120 \
		-gencode=arch=compute_120,code=compute_120
LD_FLAGS = -lcudart -L/usr/local/cuda/lib64
CPP_HEADER_FILES = inc/misc_def.h inc/host_data_structures.hpp inc/device_data_structures.h inc/file_parsing.hpp inc/simulation.h
CPP_SOURCE_FILES = src/main.cpp src/file_parsing.cpp
EXE = cuAleae

profile: ${CPP_HEADER_FILES}
	nvcc -O2 ${NVCC_INC_FLAGS} ${INC_FLAGS} -c src/simulation_kernel.cu -o simulation.o
	g++ -pg -O2 ${STD_VERSION} ${INC_FLAGS} ${CPP_SOURCE_FILES} simulation.o -o ${EXE} ${LD_FLAGS}

debug: ${CPP_HEADER_FILES}
	nvcc -G -g ${NVCC_INC_FLAGS} ${INC_FLAGS} -c src/simulation_kernel.cu -o simulation.o
	g++ -g ${STD_VERSION} ${INC_FLAGS} ${CPP_SOURCE_FILES} simulation.o -o ${EXE} ${LD_FLAGS}

build: ${CPP_HEADER_FILES}
	nvcc -O2 ${NVCC_INC_FLAGS} ${INC_FLAGS} -c src/simulation_kernel.cu -o simulation.o
	g++ -O2 ${STD_VERSION} ${INC_FLAGS} ${CPP_SOURCE_FILES} simulation.o -o ${EXE} ${LD_FLAGS}

clean:
	rm -rf $(EXE)
	rm -rf *.o