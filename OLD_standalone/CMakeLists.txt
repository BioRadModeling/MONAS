cmake_minimum_required (VERSION 3.8 FATAL_ERROR)

project (monas)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
set (CMAKE_CXX_EXTENSIONS OFF)

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(MONAS_SCRIPTS
  README.md
  )

foreach(_script ${MONAS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

message (STATUS "CMAKE_CXX_COMPILER_ID: " ${CMAKE_CXX_COMPILER_ID} )
message (STATUS "CMAKE_CXX_COMPILER_VERSION: " ${CMAKE_CXX_COMPILER_VERSION} )

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang") # using GCC || CLANG
	set (CMAKE_CXX_FLAGS_WARNINGS "-Wall -Wextra -Wpedantic")
	if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		set (CMAKE_CXX_FLAGS_WARNINGS "${CMAKE_CXX_FLAGS_WARNINGS} -Wmaybe-uninitialized")
	endif ()
	set (CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fno-trapping-math -ftree-vectorize -fno-math-errno -pthread")
	if (APPLE)
		set (CMAKE_CXX_FLAGS_BASE "")
	else ()
		set (CMAKE_CXX_FLAGS_BASE "-Wl,--disable-new-dtags")
	endif ()
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_BASE} ${CMAKE_CXX_FLAGS_WARNINGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
endif()


file (GLOB SOURCES ${PROJECT_SOURCE_DIR}/include/*.hh ${PROJECT_SOURCE_DIR}/src/*.cc)
# setup include directories for this project
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable (${PROJECT_NAME} monas.cc ${SOURCES})

# GSL libraries
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
target_link_libraries(monas ${GSL_LIBRARIES} ${GSL_CBLAS_LIBRARIES})


