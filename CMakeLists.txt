cmake_minimum_required(VERSION 3.20)
project(bnb)

set (CMAKE_CXX_STANDARD 17)

# Configuração para profilamento
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-instr-generate")

set(SOURCE_FILES
    ./src/instance-reader/Data.h
    ./src/instance-reader/Data.cpp
    ./src/hungarian/Hungarian.h
    ./src/hungarian/Hungarian.cpp
    ./src/main.cpp
    ./src/branch-and-bound/branchAndBound.h
    ./src/branch-and-bound/branchAndBound.cpp
)

add_compile_options(-O3 -Wall -Wextra -Wpedantic)
add_executable(${PROJECT_NAME} ${SOURCE_FILES})