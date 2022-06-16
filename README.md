# Getting started

Requirements
------------

* GCC 7.5
* cmake 3.16
* protobuf libraries and compiler
* OpenMPI libraries and tools
* Boost: filesystem, program_options, iostreams, date_time
* nasm

For Ubuntu 18.04 LTS the following commands installs the requirements:

`sudo apt get install --no-install-recommends g++ make libprotoc-dev protobuf-compiler libopenmpi-dev openmpi-bin libboost-filesystem-dev libboost-program-options-dev libboost-iostreams-dev libboost-date-time-dev nasm`

To get a recent cmake, download from `https://cmake.org/download/`, for example:

`wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.tar.gz`

Building
--------

`mkdir build`

`cd build`

`cmake ..`

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
* Eric Niebler. range-v3: Range library for C++14/17/20. https://github.com/ericniebler/range-v3
* Intel. ISA-L: Intelligent Storage Acceleration Library. https://github.com/intel/isa-l
* Intel. oneTBB: oneAPI Threading Building Blocks. https://github.com/oneapi-src/oneTBB
* Yatish Turakhia. UShER: Ultrafast Sample Placement on Existing Trees. https://github.com/yatisht/usher
