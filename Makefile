# Compila as versões em C++ dos métodos. Os executáveis ficam em build/, com o
# nome do arquivo de origem sem a extensão
#
# Uso: make          compila tudo
#      make clean    apaga os executáveis

CXX = g++

# -ffp-contract=off: sem fundir multiplicação e soma em uma instrução (FMA), o
# que mudaria o arredondamento em relação às versões em Python
CXXFLAGS = -O3 -march=native -std=c++17 -ffp-contract=off

EXECUTAVEIS = build/jacobi_CPU build/gauss_seidel_CPU

all: $(EXECUTAVEIS)

# As opções de compilação ficam gravadas no executável (--compilacao)
build/jacobi_CPU: jacobi/jacobi_CPU.cpp Makefile
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -DOPCOES='"$(CXXFLAGS)"' -o $@ $<

build/gauss_seidel_CPU: gauss_seidel/gauss_seidel_CPU.cpp Makefile
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -DOPCOES='"$(CXXFLAGS)"' -o $@ $<

clean:
	rm -rf build

.PHONY: all clean
