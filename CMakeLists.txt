cmake_minimum_required(VERSION 3.14)

# Some older CMake installations report the architecture of the CMake binary
# (for example x86_64 under Rosetta) instead of the native host architecture.
# Select the native Apple Silicon target before enabling the compilers.
execute_process(COMMAND uname -s
                OUTPUT_VARIABLE PTCOUNT_HOST_SYSTEM
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND uname -m
                OUTPUT_VARIABLE PTCOUNT_HOST_ARCH
                OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND sysctl -n hw.optional.arm64
                OUTPUT_VARIABLE PTCOUNT_HOST_SUPPORTS_ARM64
                OUTPUT_STRIP_TRAILING_WHITESPACE
                ERROR_QUIET)
if(PTCOUNT_HOST_SYSTEM STREQUAL "Darwin"
   AND (PTCOUNT_HOST_ARCH STREQUAL "arm64" OR PTCOUNT_HOST_SUPPORTS_ARM64 STREQUAL "1")
   AND NOT CMAKE_OSX_ARCHITECTURES)
  set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "Build architecture" FORCE)
endif()

project(PTCount LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_path(GLPK_INCLUDE_DIR glpk.h HINTS /opt/homebrew/include /usr/local/include)
find_library(GLPK_LIBRARY NAMES glpk HINTS /opt/homebrew/lib /usr/local/lib)
if(NOT GLPK_INCLUDE_DIR OR NOT GLPK_LIBRARY)
  message(FATAL_ERROR "GLPK development files are required (glpk.h and libglpk).")
endif()

set(PTCOUNT_INCLUDE_DIRS
  ${PROJECT_SOURCE_DIR}/config
  ${PROJECT_SOURCE_DIR}/graph
  ${PROJECT_SOURCE_DIR}/utility
  ${PROJECT_SOURCE_DIR}/utility/td
  ${PROJECT_SOURCE_DIR}/third_party/nauty
  ${PROJECT_SOURCE_DIR}/relation
  ${PROJECT_SOURCE_DIR}/plan
  ${PROJECT_SOURCE_DIR}/matching
  ${PROJECT_SOURCE_DIR}/counting
  ${PROJECT_SOURCE_DIR}/piso
  ${PROJECT_SOURCE_DIR}/include
  ${GLPK_INCLUDE_DIR}
)

add_library(nauty STATIC
  third_party/nauty/nauty.c
  third_party/nauty/nautil.c
  third_party/nauty/nausparse.c
  third_party/nauty/naugraph.c
  third_party/nauty/schreier.c
  third_party/nauty/naurng.c
  third_party/nauty/naugroup.c
  third_party/nauty/gtools.c
)
target_include_directories(nauty PUBLIC third_party/nauty)

add_library(td STATIC
  utility/td/td.cpp
  utility/td/fractional_cover.cpp
)
target_include_directories(td PUBLIC ${PTCOUNT_INCLUDE_DIRS})
target_link_libraries(td PUBLIC ${GLPK_LIBRARY})

add_library(ptcount_core STATIC
  utility/utils.cpp
  utility/compute_set_intersection.cpp
  graph/graph.cpp
  graph/plan_graph.cpp
  relation/candidate_space.cpp
  relation/relation.cpp
  plan/decomposition.cpp
  plan/estimator.cpp
  plan/optimizer.cpp
  plan/all_td.cpp
  plan/tree.cpp
  plan/decompose.cpp
  plan/cardinality_cache.cpp
  plan/cost_estimator.cpp
  plan/plan_cost.cpp
  matching/join.cpp
  counting/equation.cpp
  counting/execution.cpp
)
target_include_directories(ptcount_core PUBLIC ${PTCOUNT_INCLUDE_DIRS})
target_link_libraries(ptcount_core PUBLIC td nauty)

add_library(piso_components STATIC
  piso/connode.cpp
  piso/connode_plan.cpp
  piso/connode_execute.cpp
)
target_include_directories(piso_components PUBLIC ${PTCOUNT_INCLUDE_DIRS})
target_link_libraries(piso_components PUBLIC ptcount_core)

add_executable(piso src/piso.cpp)
target_link_libraries(piso PRIVATE piso_components)

add_library(piso_planner STATIC src/piso.cpp)
target_compile_definitions(piso_planner PRIVATE PISO_LIBRARY_MODE)
target_link_libraries(piso_planner PUBLIC piso_components)

add_library(tiso_runner STATIC src/tiso_runner.cpp)
target_link_libraries(tiso_runner PUBLIC ptcount_core)

add_executable(tiso src/tiso.cpp)
target_link_libraries(tiso PRIVATE tiso_runner)

add_executable(ptcount src/ptcount.cpp)
target_link_libraries(ptcount PRIVATE piso_planner tiso_runner piso_components)

set_target_properties(ptcount piso tiso PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

include(CTest)
if(BUILD_TESTING)
  find_package(Python3 COMPONENTS Interpreter)
  if(Python3_Interpreter_FOUND)
    add_test(NAME release_smoke
             COMMAND ${Python3_EXECUTABLE}
                     ${PROJECT_SOURCE_DIR}/tests/validate_release.py
                     ${CMAKE_BINARY_DIR})
  endif()
endif()
