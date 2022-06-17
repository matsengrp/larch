cmake_minimum_required(VERSION 3.16)
project (larch)

include(FetchContent)
include(ExternalProject)

find_package(Protobuf REQUIRED)
protobuf_generate_cpp(DAG_PROTO_SRCS DAG_PROTO_HDRS
    proto/dag.proto
    proto/parsimony.proto)

FetchContent_Declare(json
	URL					https://github.com/nlohmann/json/releases/download/v3.10.5/json.tar.xz
)
FetchContent_MakeAvailable(json)

ExternalProject_Add(range-v3
    GIT_REPOSITORY      https://github.com/ericniebler/range-v3.git
    GIT_TAG             0.11.0
    GIT_SHALLOW         true
    GIT_PROGRESS        true
    UPDATE_DISCONNECTED true
    CMAKE_ARGS          -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/range-v3/install -DCMAKE_BUILD_TYPE=Release -DRANGE_V3_DOCS=OFF -DRANGE_V3_TESTS=OFF -DRANGE_V3_EXAMPLES=OFF -DRANGE_V3_PERF=OFF -DRANGE_V3_HEADER_CHECKS=OFF
    PREFIX              range-v3
    INSTALL_DIR         range-v3/install
)

ExternalProject_Add(isa-l
    GIT_REPOSITORY      https://github.com/intel/isa-l.git
    GIT_TAG             v2.30.0
    GIT_SHALLOW         true
    GIT_PROGRESS        true
    UPDATE_DISCONNECTED true
    PREFIX				isa-l
    CONFIGURE_COMMAND	""
    BUILD_COMMAND		make -C ${PROJECT_BINARY_DIR}/isa-l/src/isa-l -f Makefile.unx prefix=${PROJECT_BINARY_DIR}/isa-l/install install
	INSTALL_COMMAND		""
)

ExternalProject_Add(oneTBB
    GIT_REPOSITORY      https://github.com/oneapi-src/oneTBB
	GIT_TAG             2019_U9
    GIT_SHALLOW         true
    GIT_PROGRESS        true
    UPDATE_DISCONNECTED true
    PREFIX				oneTBB
    CONFIGURE_COMMAND	""
    BUILD_COMMAND		""
	INSTALL_COMMAND		""
)

ExternalProject_Add(usher
    GIT_REPOSITORY      https://github.com/yatisht/usher
    GIT_TAG             usher-optimize-dev
    GIT_SHALLOW         true
    GIT_PROGRESS        true
    UPDATE_DISCONNECTED true
    PATCH_COMMAND       sed -i "s@oneTBB-2019_U9/include/@@g" ${PROJECT_BINARY_DIR}/usher/src/usher/src/check_samples_place/import_vcf.cpp && cat ${PROJECT_SOURCE_DIR}/add_usher_libs.cmake >> ${PROJECT_BINARY_DIR}/usher/src/usher/CMakeLists.txt
    CMAKE_CACHE_ARGS    -DCMAKE_INSTALL_PREFIX:STRING=${PROJECT_BINARY_DIR}/usher/install -DCMAKE_PREFIX_PATH:STRING=${PROJECT_BINARY_DIR}/oneTBB/src/oneTBB/cmake;${PROJECT_BINARY_DIR}/isa-l/install -DCMAKE_CXX_FLAGS:STRING=-I${PROJECT_BINARY_DIR}/isa-l/install/include
    CMAKE_ARGS			-DTBB_DIR=${PROJECT_BINARY_DIR}/oneTBB/src/oneTBB
    PREFIX              usher
    INSTALL_DIR         usher/install
)

add_dependencies(usher isa-l oneTBB)

add_library(larch_autogen OBJECT
    ${DAG_PROTO_SRCS} ${DAG_PROTO_HDRS}
)
set_target_properties(larch_autogen
    PROPERTIES CXX_CLANG_TIDY ""
)

find_package(MPI REQUIRED)
find_package(Boost COMPONENTS program_options iostreams filesystem date_time REQUIRED)

function(larch_compile_opts PRODUCT)
    target_compile_options(${PRODUCT} PUBLIC -march=native -std=c++17 -Werror -Wall -Wextra -pedantic -pedantic-errors -Wno-unused)
    if (${USE_ASAN})
        target_compile_options(${PRODUCT} PUBLIC -Og -g3 -fsanitize=address -fsanitize=undefined)
    endif()
    target_include_directories(${PRODUCT} PUBLIC include tools test)
    target_include_directories(${PRODUCT} PUBLIC ${Protobuf_INCLUDE_DIR})
    target_include_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
    target_include_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/range-v3/install/include)
    target_include_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/oneTBB/src/oneTBB/include)
    target_include_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/usher/src/usher)
    add_dependencies(${PRODUCT} range-v3 usher)
endfunction()

function(larch_link_opts PRODUCT)
	target_link_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/usher/src/usher-build/tbb_cmake_build/tbb_cmake_build_subdir_release)
    target_link_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/usher/src/usher-build)
    target_link_directories(${PRODUCT} PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/isa-l/install/lib)
    if (${USE_ASAN})
        target_link_libraries(${PRODUCT} PUBLIC -lasan -lubsan)
    endif()
    target_link_libraries(${PRODUCT} PUBLIC -lz)
    target_link_libraries(${PRODUCT} PUBLIC nlohmann_json::nlohmann_json)
    target_link_libraries(${PRODUCT} PUBLIC ${Protobuf_LIBRARIES})
    target_link_libraries(${PRODUCT} PUBLIC tbb_preview -lstdc++fs usher-sampled-lib ${Boost_LIBRARIES} ${MPI_CXX_LIBRARIES} ${MPI_CXX_LINK_FLAGS} isal)
endfunction()

function(larch_executable PRODUCT)
    add_executable(${PRODUCT} ${ARGN} $<TARGET_OBJECTS:larch_autogen>)
    larch_compile_opts(${PRODUCT})
    larch_link_opts(${PRODUCT})
    target_link_libraries(${PRODUCT} PUBLIC larch)
endfunction()

add_library(larch
    src/compact_genome.cpp
    src/dag.cpp
    src/dag_loader.cpp
    src/edge_label.cpp
    src/edge_mutations.cpp
    src/edge_storage.cpp
    src/edge.cpp
    src/leaf_set.cpp
    src/merge.cpp
    src/mutation_annotated_dag.cpp
    src/node_label.cpp
    src/node_storage.cpp
    src/node.cpp
    src/post_order_iterator.cpp
    src/pre_order_iterator.cpp)
larch_compile_opts(larch)
larch_link_opts(larch)

larch_executable(larch-test
    test/main.cpp
    test/synthetic_dags.cpp
    test/test_compact_genome.cpp
    test/test_merge.cpp
    test/test_synthetic.cpp)

larch_executable(merge
    tools/merge.cpp)

larch_executable(dag_diff
    tools/dag_diff.cpp)

larch_executable(dag2dot
    tools/dag2dot.cpp)

larch_executable(sampled
    tools/sampled.cpp)
