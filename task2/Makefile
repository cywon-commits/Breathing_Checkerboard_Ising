CXX ?= g++
FLAGS ?= -O3 -march=native -std=c++17 -fopenmp -shared -fPIC
all: libtm.so libtm_sector.so
libtm.so: src/tm.cpp
	$(CXX) $(FLAGS) $< -o $@
libtm_sector.so: src/tm_sector.cpp
	$(CXX) $(FLAGS) $< -o $@
test: all
	python3 tests/test_tm.py
