# set the version from VERSION file, and the build commit and date from git repo

message(STATUS "cmake source directory ${CMAKE_SOURCE_DIR}")

file(READ "${CMAKE_SOURCE_DIR}/VERSION" VERSION_NUMBER)
string(STRIP "${VERSION_NUMBER}" VERSION_NUMBER)

execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

execute_process(
    COMMAND git log -1 --format=%cd --date=short
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_DATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
    "tools/version.hpp.in"
    "tools/version.hpp"
    @ONLY
)
