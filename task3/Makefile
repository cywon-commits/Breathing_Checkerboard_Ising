CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -fopenmp -Wall
all: ../build/ice_pt ../build/test_ice
../build/ice_pt: src/ice_pt.cpp src/ice.hpp
	@mkdir -p ../build
	$(CXX) $(CXXFLAGS) $< -o $@
../build/test_ice: tests/test_ice.cpp src/ice.hpp
	@mkdir -p ../build
	$(CXX) $(CXXFLAGS) $< -o $@
test: ../build/test_ice
	../build/test_ice
.PHONY: all test
