# Compila as versões em C++ dos métodos. Os executáveis ficam em build/, com o
# nome do arquivo de origem sem a extensão
#
# Uso: make          compila tudo
#      make clean    apaga os executáveis

CXX = g++

# -ffp-contract=off: sem fundir multiplicação e soma em uma instrução (FMA), o
# que mudaria o arredondamento em relação às versões em Python
CXXFLAGS = -O3 -march=native -std=c++17 -ffp-contract=off

# Versões com OpenMP: as mesmas opções, mais a que ativa as diretivas
OMPFLAGS = $(CXXFLAGS) -fopenmp

# Versões em CUDA: o nvcc compila os kernels e entrega o resto ao g++, com as
# mesmas opções (-Xcompiler). --fmad=false: o mesmo que -ffp-contract=off, nos
# kernels
NVCC = /usr/local/cuda/bin/nvcc
CUDAFLAGS = -O3 -std=c++17 --fmad=false -Xcompiler -march=native,-ffp-contract=off

EXECUTAVEIS = build/jacobi_CPU build/gauss_seidel_CPU build/gauss_seidel_red_black_CPU build/successive_over_relaxation_red_black_CPU build/jacobi_CPU_openmp build/gauss_seidel_CPU_openmp build/gauss_seidel_red_black_CPU_openmp build/successive_over_relaxation_red_black_CPU_openmp build/jacobi_GPU build/gauss_seidel_GPU build/gauss_seidel_red_black_GPU build/successive_over_relaxation_red_black_GPU

all: $(EXECUTAVEIS)

# As opções de compilação ficam gravadas no executável (--compilacao)
build/jacobi_CPU: jacobi/jacobi_CPU.cpp Makefile
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -DOPCOES='"$(CXXFLAGS)"' -o $@ $<

build/gauss_seidel_CPU: gauss_seidel/gauss_seidel_CPU.cpp Makefile
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -DOPCOES='"$(CXXFLAGS)"' -o $@ $<

build/gauss_seidel_red_black_CPU: gauss_seidel_red_black/gauss_seidel_red_black_CPU.cpp Makefile
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -DOPCOES='"$(CXXFLAGS)"' -o $@ $<

build/successive_over_relaxation_red_black_CPU: successive_over_relaxation_red_black/successive_over_relaxation_red_black_CPU.cpp Makefile
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -DOPCOES='"$(CXXFLAGS)"' -o $@ $<

build/jacobi_CPU_openmp: jacobi/jacobi_CPU_openmp.cpp Makefile
	@mkdir -p build
	$(CXX) $(OMPFLAGS) -DOPCOES='"$(OMPFLAGS)"' -o $@ $<

build/gauss_seidel_CPU_openmp: gauss_seidel/gauss_seidel_CPU_openmp.cpp Makefile
	@mkdir -p build
	$(CXX) $(OMPFLAGS) -DOPCOES='"$(OMPFLAGS)"' -o $@ $<

build/gauss_seidel_red_black_CPU_openmp: gauss_seidel_red_black/gauss_seidel_red_black_CPU_openmp.cpp Makefile
	@mkdir -p build
	$(CXX) $(OMPFLAGS) -DOPCOES='"$(OMPFLAGS)"' -o $@ $<

build/successive_over_relaxation_red_black_CPU_openmp: successive_over_relaxation_red_black/successive_over_relaxation_red_black_CPU_openmp.cpp Makefile
	@mkdir -p build
	$(CXX) $(OMPFLAGS) -DOPCOES='"$(OMPFLAGS)"' -o $@ $<

build/jacobi_GPU: jacobi/jacobi_GPU.cu Makefile
	@mkdir -p build
	$(NVCC) $(CUDAFLAGS) -DOPCOES='"$(CUDAFLAGS)"' -o $@ $<

build/gauss_seidel_GPU: gauss_seidel/gauss_seidel_GPU.cu Makefile
	@mkdir -p build
	$(NVCC) $(CUDAFLAGS) -DOPCOES='"$(CUDAFLAGS)"' -o $@ $<

build/gauss_seidel_red_black_GPU: gauss_seidel_red_black/gauss_seidel_red_black_GPU.cu Makefile
	@mkdir -p build
	$(NVCC) $(CUDAFLAGS) -DOPCOES='"$(CUDAFLAGS)"' -o $@ $<

build/successive_over_relaxation_red_black_GPU: successive_over_relaxation_red_black/successive_over_relaxation_red_black_GPU.cu Makefile
	@mkdir -p build
	$(NVCC) $(CUDAFLAGS) -DOPCOES='"$(CUDAFLAGS)"' -o $@ $<

clean:
	rm -rf build

.PHONY: all clean
