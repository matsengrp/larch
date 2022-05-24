cmake_minimum_required(VERSION 3.11)
project (larch)

find_package(Protobuf REQUIRED)
protobuf_generate_cpp(DAG_PROTO_SRCS DAG_PROTO_HDRS
    proto/dag.proto
    proto/parsimony.proto)

include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.10.5/json.tar.xz)
FetchContent_MakeAvailable(json)

include(ExternalProject)
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

add_library(larch_autogen OBJECT
    ${DAG_PROTO_SRCS} ${DAG_PROTO_HDRS}
)
set_target_properties(larch_autogen
    PROPERTIES CXX_CLANG_TIDY ""
)

set(CMAKE_CXX_CLANG_TIDY clang-tidy)

function(larch_compile_opts PRODUCT)
    target_compile_options(${PRODUCT} PUBLIC -std=c++17 -Werror -Wall -Wextra -pedantic -pedantic-errors -Wno-unused)
    if (${CMAKE_BUILD_TYPE} STREQUAL "San")
        target_compile_options(${PRODUCT} PUBLIC -fsanitize=address -fsanitize=undefined)
    endif()
    target_include_directories(${PRODUCT} PUBLIC include tools test ${Protobuf_INCLUDE_DIR} ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/range-v3/install/include)
    add_dependencies(${PRODUCT} range-v3)
endfunction()

function(larch_link_opts PRODUCT)
    if (${CMAKE_BUILD_TYPE} STREQUAL "San")
        target_link_libraries(${PRODUCT} PUBLIC -lasan -lubsan)
    endif()
    target_link_libraries(${PRODUCT} PUBLIC -lz)
    target_link_libraries(${PRODUCT} PUBLIC nlohmann_json::nlohmann_json)
    target_link_libraries(${PRODUCT} PUBLIC ${Protobuf_LIBRARIES})
    target_link_libraries(${PRODUCT} PUBLIC tbb -lstdc++fs)
endfunction()

function(larch_executable PRODUCT)
    add_executable(${PRODUCT} ${ARGN} $<TARGET_OBJECTS:larch_autogen>)
    larch_compile_opts(${PRODUCT})
    larch_link_opts(${PRODUCT})
    target_link_libraries(${PRODUCT} PUBLIC larch)
endfunction()

add_library(larch
    src/history_dag.cpp
    src/history_dag_loader.cpp
    src/history_dag_node_storage.cpp
    src/zlib_stream.cpp)
larch_compile_opts(larch)
larch_link_opts(larch)

larch_executable(larch-test
    test/main.cpp
    test/test_merge.cpp)

larch_executable(merge
    tools/merge.cpp)
