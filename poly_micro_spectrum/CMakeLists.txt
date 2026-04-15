cmake_minimum_required(VERSION 3.16)
project(poly_micro_spectrum LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(poly_micro_spectrum
    src/main.cpp
    src/LookupTable.cpp
    src/LookupLibrary.cpp
    src/PhaseSpaceReader.cpp
    src/CsvWriter.cpp
    src/SpectrumAccumulator.cpp
)

target_include_directories(poly_micro_spectrum PRIVATE include)