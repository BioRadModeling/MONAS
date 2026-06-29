cmake_minimum_required(VERSION 3.16)
project(microdosimetry_with_LUTs LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(microdosimetry_with_LUTs
    src/main.cpp
    src/AmfResultParser.cpp
    src/AmfGeometryDeriver.cpp
    src/AmfRunner.cpp
    src/AmfTemplateWriter.cpp
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
    src/SpectrumMoments.cpp
)

target_include_directories(microdosimetry_with_LUTs PRIVATE include)

enable_testing()

find_package(Python3 COMPONENTS Interpreter)
if(Python3_Interpreter_FOUND)
    add_test(
        NAME amf_diagnostics
        COMMAND ${Python3_EXECUTABLE} -m unittest discover -s tests
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    )
endif()
