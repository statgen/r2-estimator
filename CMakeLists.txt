cmake_minimum_required(VERSION 3.10)
project(r2-estimator)

set(CMAKE_CXX_STANDARD 11)

if (NOT BUILD_SHARED_LIBS)
    set(CMAKE_FIND_LIBRARY_SUFFIXES ".a;${CMAKE_FIND_LIBRARY_SUFFIXES}") # Prefer libz.a when both are available
endif()

find_library(SAVVY_LIB savvy)
find_library(Z_LIB z)
find_library(ZSTD_LIB zstd)
find_library(HTS_LIB hts)

add_executable(r2-estimator src/main.cpp)
target_link_libraries(r2-estimator ${SAVVY_LIB} ${Z_LIB} ${ZSTD_LIB} ${HTS_LIB})

install(TARGETS r2-estimator RUNTIME DESTINATION bin)