# Getting started

Requirements
------------

For Ubuntu 22.04 LTS the following commands installs the requirements:

`sudo apt install --no-install-recommends git cmake make g++-10 mpi-default-dev libprotobuf-dev libboost-dev libboost-program-options-dev libboost-filesystem-dev libboost-date-time-dev protobuf-compiler pkg-config pybind11-dev libcurl4-openssl-dev python3.10-dev libhts-dev libtabixpp-dev libtabixpp0 libbz2-dev`

To get a recent cmake, download from `https://cmake.org/download/`, for example:

`wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.tar.gz`

Building
--------

`mkdir build`

`cd build`

`cmake -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER="g++-10" ..`

`make -j16`

Optionally add -DCMAKE_CXX_CLANG_TIDY="clang-tidy" to enable clang-tidy.

Optionally add -DUSE_ASAN=yes to enable asan and ubsan.

Running
-------

From the build directory:

`ln -s ../data`

`./larch-test`

Passing *nocatch* to the tests executable will allow exceptions to escape, which is useful for debugging. A gdb session can be started with `gdb --args build/larch-test nocatch`.

Third-party
-----------

* Lohmann, N. (2022). JSON for Modern C++ (Version 3.10.5) [Computer software]. https://github.com/nlohmann
* Eric Niebler. Range library for C++14/17/20. https://github.com/ericniebler/range-v3
