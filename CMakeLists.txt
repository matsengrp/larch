cmake_minimum_required(VERSION 3.16)
project(larch)

if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE yes)
set(CMAKE_EXPORT_COMPILE_COMMANDS yes)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
endif()
option(USE_USHER "Use matOptimize" yes)
option(USE_MAT_VIEW "Use MATView" yes)
option(USE_CPPTRACE "Use cpptrace" no)

if(USE_USHER)
    set(TBB_VERSION "2019_U9")
else()
    set(TBB_VERSION "v2021.11.0")
endif()

message("cmake system architecture: ${CMAKE_SYSTEM_PROCESSOR}")
message("build type: ${CMAKE_BUILD_TYPE}")
message("using usher: ${USE_USHER}")

include(FetchContent)
if(USE_USHER)
  FetchContent_Declare(oneTBB-src
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
    GIT_TAG ${TBB_VERSION}
    GIT_PROGRESS true
  )
  FetchContent_MakeAvailable(oneTBB-src)

  include(${onetbb-src_SOURCE_DIR}/cmake/TBBBuild.cmake)
  tbb_build(TBB_ROOT ${onetbb-src_SOURCE_DIR} CONFIG_DIR TBB_DIR MAKE_ARGS tbb_cpf=1 CFLAGS="-w" CXXFLAGS="-Wno-changes-meaning")

  find_package(TBB REQUIRED tbbmalloc tbbmalloc_proxy tbb_preview)
  find_package(MPI REQUIRED)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
else()
  FetchContent_Declare(oneTBB
    GIT_REPOSITORY https://github.com/oneapi-src/oneTBB
    GIT_TAG ${TBB_VERSION}
    GIT_SHALLOW true
    GIT_PROGRESS true
    UPDATE_DISCONNECTED true
    CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release -DTBB_TEST=OFF -DTBB_FIND_PACKAGE=OFF
  )
  FetchContent_MakeAvailable(oneTBB)
endif()

include(protobuf.cmake)
# find_package(Protobuf REQUIRED)

add_library(larch-autogen OBJECT
  ${DAG_PROTO_SRCS} ${DAG_PROTO_HDRS}
  ${PROTO_SRCS} ${PROTO_HDRS}
  ${DETAILED_MUTATIONS_PROTO_SRCS} ${DETAILED_MUTATIONS_PROTO_HDRS}
)
set(PROTO_OUT_DIR ${PROJECT_BINARY_DIR})
target_include_directories(larch-autogen PUBLIC ${PROTO_OUT_DIR} ${Protobuf_INCLUDE_DIRS})
target_link_libraries(larch-autogen PUBLIC ${Protobuf_ABSL_LIBS})
set_target_properties(larch-autogen
  PROPERTIES CXX_CLANG_TIDY ""
)

protobuf_generate(
  LANGUAGE cpp
  TARGET larch-autogen
  PROTOS proto/dag.proto deps/usher/parsimony.proto deps/usher/mutation_detailed.proto
)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.10.5/json.tar.xz)
FetchContent_MakeAvailable(json)

include(ExternalProject)
ExternalProject_Add(range-v3
  GIT_REPOSITORY https://github.com/ericniebler/range-v3.git
  GIT_TAG 0.12.0
  GIT_SHALLOW true
  GIT_PROGRESS true
  GIT_SUBMODULES ""
  UPDATE_DISCONNECTED true
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/range-v3/install -DCMAKE_BUILD_TYPE=Release -DRANGE_V3_DOCS=OFF -DRANGE_V3_TESTS=OFF -DRANGE_V3_EXAMPLES=OFF -DRANGE_V3_PERF=OFF -DRANGE_V3_HEADER_CHECKS=OFF
  PREFIX range-v3
  INSTALL_DIR range-v3/install
)

if(USE_CPPTRACE)
FetchContent_Declare(cpptrace
  GIT_REPOSITORY https://github.com/jeremy-rifkin/cpptrace.git
  GIT_TAG v0.8.2
  GIT_PROGRESS true)
FetchContent_MakeAvailable(cpptrace)
endif()

if(${USE_HTSLIB})
  ExternalProject_Add(htslib
    GIT_REPOSITORY https://github.com/samtools/htslib.git
    GIT_TAG 1.17
    GIT_SHALLOW true
    GIT_PROGRESS true
    UPDATE_DISCONNECTED true
    PREFIX htslib
    BUILD_IN_SOURCE true
    CONFIGURE_COMMAND autoreconf -i COMMAND ./configure --prefix=${PROJECT_BINARY_DIR}/htslib/install
    BUILD_COMMAND make -j
    INSTALL_COMMAND make install
  )
endif()

add_custom_target(generate_version_h ALL
  COMMAND cd ${CMAKE_SOURCE_DIR} && ${CMAKE_COMMAND} -DBUILD_TYPE=${CMAKE_BUILD_TYPE} -P version.cmake  
  COMMENT "Generating version.hpp with version information"
  VERBATIM
)

find_package(ZLIB)

set(STRICT_WARNINGS -Werror -fno-common -Wno-unknown-warning-option -Wno-pragmas -Wall -Wextra -pedantic -Wold-style-cast -Wshadow -Wconversion -Wcast-align -Wcast-qual -Wlogical-op -Wlogical-not-parentheses -Wredundant-decls -Wunreachable-code -Wparentheses -Wno-dangling-reference -Wno-ignored-optimization-argument -Wno-deprecated-declarations -Wno-unreachable-code-loop-increment -Wno-unneeded-internal-declaration -Wno-sign-conversion -Wno-braced-scalar-init -Wno-c++20-extensions -Wno-terminate -Wno-inaccessible-base -ftemplate-backtrace-limit=0)

if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 8.0)
  list(APPEND STRICT_WARNINGS -Wcast-align=strict)
endif()

function(larch_compile_opts PRODUCT)
  target_compile_options(${PRODUCT} PUBLIC -march=native -std=c++17 -pipe)

  if(CMAKE_BUILD_TYPE EQUAL "DEBUG")
    target_compile_options(${PRODUCT} PUBLIC -fno-omit-frame-pointer -fno-stack-protector)
  endif()

  if(${USE_USHER})
    target_compile_options(${PRODUCT} PUBLIC -DUSE_USHER)
  endif()

  if(${USE_MAT_VIEW})
    target_compile_options(${PRODUCT} PUBLIC -DUSE_MAT_VIEW)
  endif()

  if(${USE_CPPTRACE})
    target_compile_options(${PRODUCT} PUBLIC -DUSE_CPPTRACE)
  endif()

  if(${USE_ASAN})
    target_compile_options(${PRODUCT} PUBLIC -O0 -g3 -fsanitize=address,undefined -fno-sanitize-recover)
  elseif(${USE_TSAN})
    target_compile_options(${PRODUCT} PUBLIC -O0 -g3 -fsanitize=thread)
  endif()

  target_include_directories(${PRODUCT} PUBLIC include tools test)
  target_include_directories(${PRODUCT} PUBLIC ${Protobuf_INCLUDE_DIR})
  target_include_directories(${PRODUCT} PUBLIC ${PROTO_OUT_DIR})
  target_include_directories(${PRODUCT} PUBLIC ${PROTO_OUT_DIR}/deps/usher)
  target_include_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/range-v3/install/include)
  add_dependencies(${PRODUCT} range-v3)

  if(${USE_HTSLIB})
    target_include_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/htslib/install/include)
    add_dependencies(${PRODUCT} htslib)
  endif()
endfunction()

FIND_PACKAGE(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)

if(${USE_USHER})
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

function(larch_link_opts PRODUCT)
  target_link_libraries(${PRODUCT} PUBLIC larch-autogen)

  if(${USE_ASAN})
    target_link_libraries(${PRODUCT} PUBLIC -lasan -lubsan)
  elseif(${USE_TSAN})
    target_link_libraries(${PRODUCT} PUBLIC -ltsan)
  endif()

  target_link_libraries(${PRODUCT} PUBLIC -lz)
  target_link_libraries(${PRODUCT} PUBLIC nlohmann_json::nlohmann_json)
  target_link_libraries(${PRODUCT} PUBLIC ${Protobuf_LIBRARIES})
  target_link_libraries(${PRODUCT} PUBLIC ${Boost_LIBRARIES})

  if(${USE_USHER})
    target_link_libraries(${PRODUCT} PUBLIC ${TBB_IMPORTED_TARGETS} -lstdc++fs)
  else()
    target_link_libraries(${PRODUCT} PUBLIC TBB::tbb)
  endif()

  if(${USE_HTSLIB})
    target_link_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/htslib/install/lib)
    target_link_libraries(${PRODUCT} PUBLIC -lhts)
  endif()

  if(${USE_CPPTRACE})
    target_link_libraries(${PRODUCT} PUBLIC cpptrace::cpptrace)
  endif()
endfunction()

function(larch_executable PRODUCT)
  add_executable(${PRODUCT} ${ARGN})
  larch_compile_opts(${PRODUCT})
  target_compile_options(${PRODUCT} PRIVATE ${STRICT_WARNINGS})
  larch_link_opts(${PRODUCT})
  target_link_libraries(${PRODUCT} PUBLIC larch)
  add_dependencies(${PRODUCT} generate_version_h)
endfunction()

add_library(larch
  src/dag_loader.cpp
)
larch_compile_opts(larch)
target_compile_options(larch PRIVATE ${STRICT_WARNINGS})

if(${USE_ASAN})
else()
  target_compile_options(larch PRIVATE -Wstack-usage=2048)
endif()

larch_link_opts(larch)

if(${USE_USHER})
  add_library(usher
    ${MATOPTIMIZE_SRC}
    src/usher_globals.cpp)
  target_compile_options(usher PRIVATE -march=native -std=c++17 -w -fno-omit-frame-pointer -DUSE_USHER)

  if(${USE_ASAN})
    target_compile_options(usher PUBLIC -O0 -g3 -fsanitize=address,undefined -fno-sanitize-recover)
  elseif(${USE_TSAN})
    target_compile_options(usher PUBLIC -O0 -g3 -fsanitize=thread)
  endif()

  target_include_directories(usher PUBLIC include)
  target_include_directories(usher PUBLIC ${Protobuf_INCLUDE_DIR})
  target_include_directories(usher PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/range-v3/install/include)

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

if(${USE_USHER})
  target_link_libraries(larch-test PUBLIC larch-usher-glue usher)
  add_dependencies(larch-test larch-usher-glue)
  set_target_properties(larch-test
    PROPERTIES CXX_CLANG_TIDY ""
  )
endif()

larch_executable(larch-dagutil
  tools/larch-dagutil.cpp)

larch_executable(larch-dag2dot
  tools/larch-dag2dot.cpp)

if(${USE_USHER})
  larch_executable(larch-usher
    tools/larch-usher.cpp
  )
  target_compile_options(larch-usher PRIVATE ${STRICT_WARNINGS})
  target_link_libraries(larch-usher PRIVATE larch-usher-glue)
  add_dependencies(larch-usher larch-usher-glue)
endif()
