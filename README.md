# Getting started

Requirements
------------

* GCC 7.5
* cmake 3.16

For Ubuntu 18.04 LTS the following commands installs the requirements:

`sudo apt install --no-install-recommends git cmake make g++ mpi-default-dev libprotobuf-dev libboost-dev libboost-program-options-dev libboost-filesystem-dev libboost-iostreams-dev libboost-date-time-dev protobuf-compiler automake autoconf libtool nasm`

To get a recent cmake, download from `https://cmake.org/download/`, for example:

`wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.tar.gz`

To setup a conda environment capable of building Larch, use the environment
file provided:

`conda env create -f environment.yml`


Building
--------

`mkdir build`

`cd build`

`cmake -DCMAKE_BUILD_TYPE=Debug ..`

`make -j16`

Optionally add -DCMAKE_CXX_CLANG_TIDY="clang-tidy" to enable clang-tidy.

Optionally add -DUSE_ASAN=yes to enable asan and ubsan.

Running
-------

From the build directory:

`ln -s ../data`

`./larch-test`

Passing *nocatch* to the tests executable will allow exceptions to escape, which is useful for debugging. A gdb session can be started with `gdb --args build/larch-test nocatch`.

*--list* produces a list of all available tests, along with an ID number.
*--range* runs tests by ID in [begin, end] range arguments.

Third-party
-----------

* Lohmann, N. (2022). JSON for Modern C++ (Version 3.10.5) [Computer software]. https://github.com/nlohmann
* Eric Niebler. Range library for C++14/17/20. https://github.com/ericniebler/range-v3
