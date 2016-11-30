cmake_minimum_required(VERSION 3.6)
project(SuperComputing)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp ProcessConfigurator.h FiniteDifference.cpp FiniteDifference.h Condition.h Var6Cond.h ResultWritter.cpp ResultWritter.h ProcessConfigurator.cpp)
add_executable(SuperComputing ${SOURCE_FILES})