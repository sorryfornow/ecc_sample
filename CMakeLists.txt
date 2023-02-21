cmake_minimum_required(VERSION 3.23)
project(ecc_sample)

set(CMAKE_CXX_STANDARD 23)

add_executable(ecc_sample playground.cpp DYGraph_Unified.hpp DYGraph_Base.hpp definition.h)
