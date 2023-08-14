include(FetchContent)

set(Protobuf_USE_STATIC_LIBS ON)
cmake_policy(SET CMP0026 OLD)

FetchContent_Declare(
  protocolbuffers_protobuf
  GIT_REPOSITORY https://github.com/protocolbuffers/protobuf.git
  GIT_TAG v23.4
  OVERRIDE_FIND_PACKAGE
)
FetchContent_MakeAvailable(protocolbuffers_protobuf)

set(Protobuf_ROOT ${protocolbuffers_protobuf_SOURCE_DIR})

message(STATUS "Setting up protobuf ...")
execute_process(
  COMMAND
    ${CMAKE_COMMAND} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -D protobuf_BUILD_TESTS=OFF -D protobuf_BUILD_PROTOC_BINARIES=ON -D CMAKE_POSITION_INDEPENDENT_CODE=ON -G "${CMAKE_GENERATOR}" .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${Protobuf_ROOT})
if(result)
  message(FATAL_ERROR "Failed to download protobuf (${result})!")
endif()

message(STATUS "Building protobuf ...")
execute_process(
  COMMAND ${CMAKE_COMMAND} --build .
  RESULT_VARIABLE result
  WORKING_DIRECTORY ${Protobuf_ROOT})
if(result)
  message(FATAL_ERROR "Failed to build protobuf (${result})!")
endif()

message(STATUS "Installing protobuf ...")
execute_process(
    COMMAND ${CMAKE_COMMAND} --install . --config ${CMAKE_BUILD_TYPE} --prefix ${Protobuf_ROOT}/install
RESULT_VARIABLE result
WORKING_DIRECTORY ${Protobuf_ROOT})
if(result)
    message(FATAL_ERROR "Failed to build protobuf (${result})!")
endif()

if (NOT EXISTS ${Protobuf_ROOT}/install/lib64)
  file (CREATE_LINK ${Protobuf_ROOT}/install/lib ${Protobuf_ROOT}/install/lib64 SYMBOLIC)
endif()

find_package(Protobuf CONFIG REQUIRED HINTS ${Protobuf_ROOT}/install/lib64/cmake)

set(Protobuf_PROTOC_EXECUTABLE ${Protobuf_ROOT}/install/bin/protoc)

include(${Protobuf_ROOT}/install/lib64/cmake/protobuf/protobuf-config.cmake)
include(${Protobuf_ROOT}/install/lib64/cmake/protobuf/protobuf-module.cmake)
include(${Protobuf_ROOT}/install/lib64/cmake/protobuf/protobuf-options.cmake)
include(${Protobuf_ROOT}/install/lib64/cmake/protobuf/protobuf-targets.cmake)

if(Protobuf_FOUND)
  message(STATUS "Protobuf version : ${Protobuf_VERSION}")
  message(STATUS "Protobuf include path : ${Protobuf_INCLUDE_DIRS}")
  message(STATUS "Protobuf libraries : ${Protobuf_LIBRARIES}")
  message(STATUS "Protobuf compiler libraries : ${Protobuf_PROTOC_LIBRARIES}")
  message(STATUS "Protobuf lite libraries : ${Protobuf_LITE_LIBRARIES}")
  message(STATUS "Protobuf protoc : ${Protobuf_PROTOC_EXECUTABLE}")
else()
  message(
    FATAL_ERROR
      "Protobuf package not found -> specify search path via Protobuf_ROOT variable"
  )
endif()

set(Protobuf_ABSL_LIBS absl::absl_check absl::absl_log absl::algorithm absl::base absl::bind_front absl::bits absl::btree absl::cleanup absl::cord absl::core_headers absl::debugging absl::die_if_null absl::dynamic_annotations absl::flags absl::flat_hash_map absl::flat_hash_set absl::function_ref absl::hash absl::layout absl::log_initialize absl::log_severity absl::memory absl::node_hash_map absl::node_hash_set absl::optional absl::span absl::status absl::statusor absl::strings absl::synchronization absl::time absl::type_traits absl::utility absl::variant utf8_range utf8_validity)
