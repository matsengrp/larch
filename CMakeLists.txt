cmake_minimum_required(VERSION 3.16)
project(larch)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  message(WARNING "ERROR: CMake build directory cannot be the same as root source directory.")
endif()

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib/static)
set(CMAKE_DEPS_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/deps)

set(CMAKE_SHARED_LIBRARY_PREFIX "liblarch_")
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib;${CMAKE_INSTALL_PREFIX}/lib/larch")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(USE_USHER "Use matOptimize" ON)
option(USE_MAT_VIEW "Use MATView" ON)
option(USE_CPPTRACE "Use cpptrace" OFF)
option(USE_ASAN "Use asan" OFF)
option(USE_TSAN "Use tsan" OFF)
option(USE_HTSLIB "Use htslib" OFF)

# set default cmake build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
endif()

set(TBB_VERSION "v2022.1.0")

message(STATUS "cmake system architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "cmake build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "using usher: ${USE_USHER}")

# External static libraries: protobuf/larch-autogen, range-v3, nlohmann_json, htslib, oneTBB, boost, zlib
include(FetchContent)
include(ExternalProject)

function(lib_install PRODUCT)
  get_target_property(LIB_PATH ${PRODUCT} IMPORTED_LOCATION_RELEASE)
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    get_target_property(LIB_PATH ${PRODUCT} IMPORTED_LOCATION_DEBUG)
  endif()
  get_filename_component(LIB_NAME ${LIB_PATH} NAME)
  install(FILES "${LIB_PATH}"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/larch"
  )
endfunction()

# # oneTBB
FetchContent_Declare(oneTBB
  GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
  GIT_TAG ${TBB_VERSION}
  GIT_SHALLOW true
  GIT_PROGRESS true
  UPDATE_DISCONNECTED true
  CMAKE_ARGS -DTBB_TEST:BOOL=OFF
)
set(TBB_TEST OFF CACHE INTERNAL "" FORCE)
FetchContent_MakeAvailable(oneTBB)
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# # protobuf
find_package(Protobuf REQUIRED)

add_library(larch-autogen OBJECT
  ${DAG_PROTO_SRCS} ${DAG_PROTO_HDRS}
  ${PROTO_SRCS} ${PROTO_HDRS}
  ${DETAILED_MUTATIONS_PROTO_SRCS} ${DETAILED_MUTATIONS_PROTO_HDRS}
)
set(PROTO_OUT_DIR ${PROJECT_BINARY_DIR})

# # larch-autogen library
target_include_directories(larch-autogen PUBLIC ${PROTO_OUT_DIR} ${Protobuf_INCLUDE_DIRS})
target_link_libraries(larch-autogen PUBLIC ${Protobuf_LIBRARIES})
set_target_properties(larch-autogen
  PROPERTIES CXX_CLANG_TIDY ""
)
protobuf_generate(
  LANGUAGE cpp
  TARGET larch-autogen
  PROTOS proto/dag.proto deps/usher/parsimony.proto deps/usher/mutation_detailed.proto
)

# # nlohmann_json
FetchContent_Declare(json
  URL https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz
  BINARY_DIR ${CMAKE_BINARY_DIR}/deps/nlohmann_json
)
FetchContent_MakeAvailable(json)

# # range-v3
ExternalProject_Add(range-v3
  GIT_REPOSITORY https://github.com/ericniebler/range-v3.git
  GIT_TAG 0.12.0
  GIT_SHALLOW true
  GIT_PROGRESS true
  GIT_SUBMODULES ""
  UPDATE_DISCONNECTED true
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/deps/range-v3/install -DCMAKE_BUILD_TYPE=Release -DRANGE_V3_DOCS=OFF -DRANGE_V3_TESTS=OFF -DRANGE_V3_EXAMPLES=OFF -DRANGE_V3_PERF=OFF -DRANGE_V3_HEADER_CHECKS=OFF
  PREFIX deps/range-v3
  INSTALL_DIR deps/range-v3/install
)

if(USE_CPPTRACE)
FetchContent_Declare(cpptrace
  GIT_REPOSITORY https://github.com/jeremy-rifkin/cpptrace.git
  GIT_TAG v0.8.2
  GIT_PROGRESS true)
FetchContent_MakeAvailable(cpptrace)
endif()

# # htslib
if(USE_HTSLIB)
  ExternalProject_Add(htslib
    GIT_REPOSITORY https://github.com/samtools/htslib.git
    GIT_TAG 1.17
    GIT_SHALLOW true
    GIT_PROGRESS true
    UPDATE_DISCONNECTED true
    PREFIX deps/htslib
    BUILD_IN_SOURCE true
    CONFIGURE_COMMAND autoreconf -i COMMAND ./configure --prefix=${PROJECT_BINARY_DIR}/deps/htslib/install
    BUILD_COMMAND make -j
    INSTALL_COMMAND make install
  )
endif()

# # zlib
find_package(ZLIB REQUIRED)

# # boost
find_package(Boost REQUIRED COMPONENTS program_options filesystem date_time)

# Internal libraries and binaries

# # set build version
add_custom_target(generate_version_h ALL
  COMMAND cd ${CMAKE_SOURCE_DIR} && ${CMAKE_COMMAND} -DBUILD_TYPE=${CMAKE_BUILD_TYPE} -P version.cmake
  COMMENT "Generating version.hpp with version information"
  VERBATIM
)

# # set compiler warnings
set(STRICT_WARNINGS -Werror -fno-common -Wno-unknown-warning-option -Wno-pragmas -Wall -Wextra -pedantic -Wold-style-cast -Wshadow -Wconversion -Wcast-align -Wcast-qual -Wlogical-op -Wlogical-not-parentheses -Wredundant-decls -Wunreachable-code -Wparentheses -Wno-dangling-reference -Wno-ignored-optimization-argument)

# # set compiler version
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.0)
  list(APPEND STRICT_WARNINGS -Wcast-align=strict)
endif()

# set larch compiler options
function(larch_compile_opts PRODUCT)
  target_compile_options(${PRODUCT} PUBLIC -march=native -std=c++17 -pipe)

  if(CMAKE_BUILD_TYPE EQUAL "DEBUG")
    target_compile_options(${PRODUCT} PUBLIC -fno-omit-frame-pointer -fno-stack-protector)
  endif()

  if(USE_USHER)
    target_compile_options(${PRODUCT} PUBLIC -DUSE_USHER)
  endif()

  if(USE_MAT_VIEW)
    target_compile_options(${PRODUCT} PUBLIC -DUSE_MAT_VIEW)
  endif()

  if(USE_CPPTRACE)
    target_compile_options(${PRODUCT} PUBLIC -DUSE_CPPTRACE)
  endif()

  if(USE_ASAN)
    target_compile_options(${PRODUCT} PUBLIC -O0 -g3 -fsanitize=address -fno-sanitize-recover)
  elseif(USE_TSAN)
    target_compile_options(${PRODUCT} PUBLIC -O0 -g3 -fsanitize=thread)
  endif()

  target_include_directories(${PRODUCT} PUBLIC include tools test)
  target_include_directories(${PRODUCT} PUBLIC ${Protobuf_INCLUDE_DIR})
  target_include_directories(${PRODUCT} PUBLIC ${PROTO_OUT_DIR})
  target_include_directories(${PRODUCT} PUBLIC ${PROTO_OUT_DIR}/deps/usher)
  target_include_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/deps/range-v3/install/include)
  add_dependencies(${PRODUCT} range-v3)

  if(USE_HTSLIB)
    target_include_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/deps/htslib/install/include)
    add_dependencies(${PRODUCT} htslib)
  endif()
endfunction()

if(USE_USHER)
  include_directories(deps/usher/)
  file(GLOB MATOPTIMIZE_SRC
    "deps/usher/src/matOptimize/apply_move/*.cpp"
    "deps/usher/src/matOptimize/Profitable_Moves_Enumerators/*.cpp"
    "deps/usher/src/matOptimize/Fitch_Sankoff.cpp"
    "deps/usher/src/matOptimize/mutation_annotated_tree.cpp"
    "deps/usher/src/matOptimize/mutation_annotated_tree_nuc_util.cpp"
    "deps/usher/src/matOptimize/mutation_annotated_tree_node.cpp"
    "deps/usher/src/matOptimize/mutation_annotated_tree_load_store.cpp"
    "deps/usher/src/matOptimize/mutation_annotated_tree_nuc_utils.cpp"
    "deps/usher/src/matOptimize/optimize_inner_loop.cpp"
    "deps/usher/src/matOptimize/Mutation_Collection.cpp"
    "deps/usher/src/matOptimize/condense.cpp"
    "deps/usher/src/matOptimize/reassign_states.cpp"
    "deps/usher/src/matOptimize/detailed_mutations_load.cpp"
    "deps/usher/src/matOptimize/detailed_mutations_store.cpp"
    "deps/usher/src/matOptimize/optimize_tree.cpp"
    "deps/usher/src/matOptimize/priority_conflict_resolver.cpp"
    "deps/usher/src/matOptimize/check_samples.cpp"
  )
endif()

# link libraries to executable
function(larch_link_opts PRODUCT)
  # Workaround for protobuf
  target_link_options(${PRODUCT} PUBLIC -Wl,--copy-dt-needed-entries)

  target_link_libraries(${PRODUCT} PUBLIC larch-autogen)

  if(USE_ASAN)
    target_link_libraries(${PRODUCT} PUBLIC -lasan)
  elseif(USE_TSAN)
    target_link_libraries(${PRODUCT} PUBLIC -ltsan)
  endif()

  target_link_libraries(${PRODUCT} PUBLIC -lz)
  target_link_libraries(${PRODUCT} PUBLIC nlohmann_json::nlohmann_json)
  target_link_libraries(${PRODUCT} PUBLIC ${Protobuf_LIBRARIES})
  target_link_libraries(${PRODUCT} PUBLIC ${Boost_LIBRARIES})

  target_link_libraries(${PRODUCT} PUBLIC TBB::tbb)

  if(USE_HTSLIB)
    target_link_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/htslib/install)
    target_link_libraries(${PRODUCT} PUBLIC -lhts)
  endif()
endfunction()

# build larch executable
function(larch_executable PRODUCT)
  add_executable(${PRODUCT} ${ARGN})
  larch_compile_opts(${PRODUCT})
  target_compile_options(${PRODUCT} PRIVATE ${STRICT_WARNINGS})
  larch_link_opts(${PRODUCT})
  target_link_libraries(${PRODUCT} PUBLIC larch)
  add_dependencies(${PRODUCT} generate_version_h)
endfunction()

function(larch_install PRODUCT)
  install(TARGETS ${PRODUCT}
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
  )
endfunction()

# Internal libraries: larch, usher

# # larch library
add_library(larch
  src/dag_loader.cpp
)
larch_compile_opts(larch)
target_compile_options(larch PRIVATE ${STRICT_WARNINGS})

if(USE_ASAN)
else()
  target_compile_options(larch PRIVATE -Wstack-usage=2048)
endif()

larch_link_opts(larch)

# # usher library
if(USE_USHER)
  add_library(usher
    ${MATOPTIMIZE_SRC}
    src/usher_globals.cpp
  )
  target_compile_options(usher PRIVATE -march=native -std=c++17 -w -fno-omit-frame-pointer -DUSE_USHER)

  if(${USE_ASAN})
    target_compile_options(usher PUBLIC -O0 -g3 -fsanitize=address -fno-sanitize-recover)
  elseif(${USE_TSAN})
    target_compile_options(usher PUBLIC -O0 -g3 -fsanitize=thread)
  endif()

  target_include_directories(usher PUBLIC include)
  target_include_directories(usher PUBLIC ${Protobuf_INCLUDE_DIR})
  target_include_directories(usher PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/deps/range-v3/install/include)

  add_dependencies(usher range-v3)

  target_include_directories(usher PUBLIC ${PROTO_OUT_DIR}/deps/usher)
  larch_link_opts(usher)
  target_link_libraries(usher PRIVATE ${Boost_LIBRARIES} ZLIB::ZLIB ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS})

  set_target_properties(usher
    PROPERTIES CXX_CLANG_TIDY ""
  )

  add_library(larch-usher-glue
    src/usher_glue.cpp)
  larch_compile_opts(larch-usher-glue)
  target_compile_options(larch-usher-glue PRIVATE ${STRICT_WARNINGS})
  larch_link_opts(larch-usher-glue)
  target_link_libraries(larch-usher-glue PUBLIC usher)
  add_dependencies(larch-usher-glue larch usher)
  set_target_properties(larch-usher-glue
    PROPERTIES CXX_CLANG_TIDY ""
  )
endif()

# Executables: larch-usher, larch-dagutil, larch-dag2dot, larch-test

# # larch-test
larch_executable(larch-test
  test/main.cpp
  test/test_ambiguous_vcf.cpp
  test/test_ambiguities.cpp
  test/test_compact_genome.cpp
  test/test_count_trees.cpp
  test/test_dag_completion.cpp
  test/test_dag_trimming.cpp
  test/test_fileio_dagbin.cpp
  test/test_fileio_protobuf.cpp
  test/test_larch_usher.cpp
  test/test_lca.cpp
  test/test_map.cpp
  test/test_mat_conversion.cpp
  test/test_mat_view.cpp
  test/test_matOptimize.cpp
  test/test_merge.cpp
  test/test_overlay.cpp
  test/test_rf_distance.cpp
  test/test_sample_tree.cpp
  test/test_sample_best_tree.cpp
  test/test_spr_after_move.cpp
  test/test_spr.cpp
  test/test_subtree_weight.cpp
  test/test_weight_accum.cpp
  test/test_weight_counter.cpp
  test/test_write_parsimony_protobuf.cpp
)
target_compile_options(larch-test PRIVATE ${STRICT_WARNINGS})

if(USE_USHER)
  target_link_libraries(larch-test PUBLIC larch-usher-glue usher)
  add_dependencies(larch-test larch-usher-glue)
  set_target_properties(larch-test
    PROPERTIES CXX_CLANG_TIDY ""
  )
endif()

# # larch-dagutil
larch_executable(larch-dagutil
  tools/larch-dagutil.cpp
)
larch_install(larch-dagutil)

# # larch-dag2dot
larch_executable(larch-dag2dot
  tools/larch-dag2dot.cpp
)
larch_install(larch-dag2dot)

# # larch-usher
if(USE_USHER)
  larch_executable(larch-usher
    tools/larch-usher.cpp
  )
  target_compile_options(larch-usher PRIVATE ${STRICT_WARNINGS})
  target_link_libraries(larch-usher PRIVATE larch-usher-glue)
  add_dependencies(larch-usher larch-usher-glue)
  larch_install(larch-usher)
endif()

