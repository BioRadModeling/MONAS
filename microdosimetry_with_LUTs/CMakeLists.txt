cmake_minimum_required(VERSION 3.16)
project(microdosimetry_with_LUTs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(microdosimetry_with_LUTs
    src/main.cpp
    src/InaniwaCalculator.cpp
    src/InaniwaLookup.cpp
    src/LetAccumulator.cpp
    src/LetCalculator.cpp
    src/LetLookup.cpp
    src/MaginiCalculator.cpp
    src/MaginiLookup.cpp
    src/LookupTable.cpp
    src/LookupLibrary.cpp
    src/PhaseSpaceReader.cpp
    src/CsvWriter.cpp
    src/SpectrumAccumulator.cpp
)

target_include_directories(microdosimetry_with_LUTs PRIVATE include)
