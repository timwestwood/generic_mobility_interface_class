CILIA_CPP = cilia_sim_main.cpp matrix.cpp quaternion.cpp segment.cpp filament.cpp broyden_solver.cpp rigid_body.cpp swimmer.cpp mobility_solver.cpp
CILIA_CUDA = cuda_functions.cu rpy_mobility_solver.cu stokesdrag_mobility_solver.cu seeding.cu

FLOW_FIELD_CPP = flow_field_main.cpp matrix.cpp quaternion.cpp
FLOW_FIELD_CUDA = flow_field_evaluator.cu


cilia_clean:
	-rm cilia
	-rm cilia.exe
	-rm cilia.exp
	-rm cilia.lib

# cilia_data_clean:
# 	-rm 

flow_field_clean:
	-rm flow_field
	-rm flow_field.exe
	-rm flow_field.exp
	-rm flow_field.lib

#
# The following work on the machines at Imperial:
#

# Basic expression works fine on nvidia2.
cilia_nvidia2:
	nvcc -g $(CILIA_CPP) $(CILIA_CUDA) -o cilia -O3 -llapack -lblas -std=c++11 -w

# Compiles but just outputs zeros if you don't specify the 3.5 architecture of nvidia3's Tesla K40 GPUs.
cilia_nvidia3:
	nvcc -g $(CILIA_CPP) $(CILIA_CUDA) -o cilia -O3 -llapack -lblas -arch=sm_35 -std=c++11 -w

flow_field_nvidia3:
	nvcc -g $(FLOW_FIELD_CPP) $(FLOW_FIELD_CUDA) -o flow_field -O3 -llapack -lblas -arch=sm_35 -std=c++11 -w

# Basic expression works fine on nvidia4.
cilia_nvidia4:
	nvcc -g $(CILIA_CPP) $(CILIA_CUDA) -o cilia -O3 -llapack -lopenblas -std=c++11 -w

flow_field_nvidia4:
	nvcc -g $(FLOW_FIELD_CPP) $(FLOW_FIELD_CUDA) -o flow_field -O3 -llapack -lopenblas -std=c++11 -w

# Basic expression works fine on HPC cluster once we've loaded the cuda module.
cilia_ic_hpc:
	module load cuda/11.4.2 && \
	nvcc -g $(CILIA_CPP) $(CILIA_CUDA) -o cilia -O3 -lopenblas -std=c++11 -w \
	-L/usr/lib64 -l:liblapack.so.3.8 # There doesn't seem to be the right symbolic link on these machines to avoid having to give the version number explicitly. Call "locate liblapack.so" to check for up-to-date versions and their locations if this fails in the future.

flow_field_ic_hpc:
	module load cuda/11.4.2 && \
	nvcc -g $(FLOW_FIELD_CPP) $(FLOW_FIELD_CUDA) -o flow_field -O3 -lopenblas -std=c++11 -w \
	-L/usr/lib64 -l:liblapack.so.3.8

#
# The following work on my personal PC. Compiling CUDA on Windows requires using the Microsoft
# Visual C++ compiler cl.exe and as such the paths etc. follow the Windows, rather than unix, format.
# blas_win64_MT.dll and lapack_win64_MT.dll have been copied into Cygwin's /bin/
#

debug_cilia_pc:
	nvcc -g $(CILIA_CPP) $(CILIA_CUDA) -o cilia -I "D:\cygwin64\lib\include" \
	-ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Tools\MSVC\14.21.27702\bin\Hostx64\x64" \
	"D:\cygwin64\lib\lib_win64\blas_win64_MT.lib" "D:\cygwin64\lib\lib_win64\lapack_win64_MT.lib" -lineinfo

cilia_pc:
	nvcc -g $(CILIA_CPP) $(CILIA_CUDA) -o cilia -I "D:\cygwin64\lib\include" \
	-ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Tools\MSVC\14.21.27702\bin\Hostx64\x64" \
	"D:\cygwin64\lib\lib_win64\blas_win64_MT.lib" "D:\cygwin64\lib\lib_win64\lapack_win64_MT.lib" -O3 -lineinfo -w

flow_field_pc:
	nvcc -g $(FLOW_FIELD_CPP) $(FLOW_FIELD_CUDA) -o flow_field -O3 \
	-ccbin "C:\Program Files (x86)\Microsoft Visual Studio\2019\BuildTools\VC\Tools\MSVC\14.21.27702\bin\Hostx64\x64" \
	"D:\cygwin64\lib\lib_win64\blas_win64_MT.lib" "D:\cygwin64\lib\lib_win64\lapack_win64_MT.lib"
