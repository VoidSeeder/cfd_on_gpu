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

EXECUTAVEIS = build/jacobi_CPU build/gauss_seidel_CPU build/gauss_seidel_red_black_CPU build/successive_over_relaxation_red_black_CPU build/jacobi_CPU_openmp

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

clean:
	rm -rf build

.PHONY: all clean
