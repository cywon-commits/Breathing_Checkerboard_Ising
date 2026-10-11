CXX ?= g++
CXXFLAGS ?= -O3 -march=native -std=c++17 -fopenmp -Wall
all: build/bcising build/test_main
build/bcising: src/bcising.cpp src/model.hpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@
build/test_main: tests/test_main.cpp src/model.hpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@
test: build/test_main
	./build/test_main
clean:
	rm -rf build
.PHONY: all test clean
