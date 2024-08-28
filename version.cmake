# set the version from VERSION file, and the build commit and date from git repo

# get git version tag
execute_process(
    COMMAND git describe --tags --abbrev=0
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_TAG_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)

# get version from git tag if available, else fall back to VERSION file
if(GIT_TAG_VERSION)
    set(VERSION_NUMBER ${GIT_TAG_VERSION})
    message("using version from git tag: ${VERSION_NUMBER}")
else()
    file(READ "${CMAKE_SOURCE_DIR}/VERSION" VERSION_NUMBER)
    string(STRIP "${VERSION_NUMBER}" VERSION_NUMBER)
    set(VERSION_NUMBER "v${VERSION_NUMBER}")
    message("using version from VERSION file: ${VERSION_NUMBER}")
endif()

# get git hash
execute_process(
    COMMAND git rev-parse --short HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# get git date
execute_process(
    COMMAND git log -1 --format=%cd --date=short
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_DATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# output version details to file
configure_file(
    "tools/version.hpp.in"
    "tools/version.hpp"
    @ONLY
)
