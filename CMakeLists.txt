cmake_minimum_required(VERSION 3.6)
project(SuperComputing)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp ProcessConfigurator.h FiniteDifference.cpp FiniteDifference.h)
add_executable(SuperComputing ${SOURCE_FILES})