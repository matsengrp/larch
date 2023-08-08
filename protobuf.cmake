include(FetchContent)

set(Protobuf_USE_STATIC_LIBS ON)

FetchContent_Declare(
  protocolbuffers_protobuf
  GIT_REPOSITORY https://github.com/protocolbuffers/protobuf.git
  GIT_TAG v23.4
  OVERRIDE_FIND_PACKAGE
)
# set(protobuf_MODULE_COMPATIBLE ON)
FetchContent_MakeAvailable(protocolbuffers_protobuf)

set(Protobuf_ROOT ${protocolbuffers_protobuf_SOURCE_DIR})
set(Protobuf_DIR ${protocolbuffers_protobuf_BINARY_DIR}/cmake/protobuf)

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
if(WIN32)
    execute_process(
        COMMAND ${CMAKE_COMMAND} --build . --target install --config ${CMAKE_BUILD_TYPE}
    RESULT_VARIABLE result
    WORKING_DIRECTORY ${Protobuf_ROOT})
    if(result)
        message(FATAL_ERROR "Failed to build protobuf (${result})!")
    endif()
endif()

find_package(Protobuf CONFIG REQUIRED HINTS ${Protobuf_DIR})

# file(TOUCH ${Protobuf_DIR}/protobuf-targets.cmake)
# include(${Protobuf_DIR}/protobuf-config.cmake)
# include(${Protobuf_DIR}/protobuf-module.cmake)
# include(${Protobuf_DIR}/protobuf-options.cmake)
# include(${Protobuf_DIR}/protobuf-targets.cmake)

if(Protobuf_FOUND)
  message(STATUS "Protobuf version : ${Protobuf_VERSION}")
  message(STATUS "Protobuf include path : ${Protobuf_INCLUDE_DIRS}")
  message(STATUS "Protobuf libraries : ${Protobuf_LIBRARIES}")
  message(STATUS "Protobuf compiler libraries : ${Protobuf_PROTOC_LIBRARIES}")
  message(STATUS "Protobuf lite libraries : ${Protobuf_LITE_LIBRARIES}")
  message(STATUS "Protobuf protoc : ${Protobuf_PROTOC_EXECUTABLE}")
else()
  message(
    WARNING
      "Protobuf package not found -> specify search path via Protobuf_ROOT variable"
  )
endif()
