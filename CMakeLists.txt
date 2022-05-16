cmake_minimum_required(VERSION 3.11)
project (larch)

include(FetchContent)
FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.10.5/json.tar.xz)
FetchContent_MakeAvailable(json)

find_package(Protobuf REQUIRED)
protobuf_generate_cpp(DAG_PROTO_SRCS DAG_PROTO_HDRS proto/dag.proto)

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

add_executable(larch
    ${DAG_PROTO_SRCS} ${DAG_PROTO_HDRS}
    src/history_dag.cpp
    src/history_dag_loader.cpp
    src/history_dag_node_storage.cpp
    src/zlib_stream.cpp
    test/main.cpp
    test/test_merge.cpp)

add_dependencies(larch range-v3)
target_compile_options(larch PUBLIC -std=c++17 -Werror -Wall -Wextra -pedantic -pedantic-errors -Wno-unused)
target_include_directories(larch PUBLIC include test ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_BINARY_DIR}/range-v3/install/include)
target_link_libraries(larch PUBLIC -lz)
target_link_libraries(larch PUBLIC nlohmann_json::nlohmann_json)
target_link_libraries(larch PUBLIC ${Protobuf_LIBRARIES})
