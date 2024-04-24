# Getting started

Requirements
------------

* GCC 7.5
* cmake 3.16

For Ubuntu 18.04 LTS the following commands installs the requirements:

`sudo apt install --no-install-recommends git cmake make g++ mpi-default-dev libprotobuf-dev libboost-dev libboost-program-options-dev libboost-filesystem-dev libboost-iostreams-dev libboost-date-time-dev protobuf-compiler automake autoconf libtool nasm`

To get a recent cmake, download from `https://cmake.org/download/`, for example:

`wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.tar.gz`

Build Environments
------------------

* singularity 3.5.3
* conda 22.9.0

Larch can be built utilizing a Singularity container or a Conda environment.

To build Singularity image, use the definition provided:

```shell
singularity build larch-singularity.sif larch-singularity.def
singularity shell larch-singularity.sif --net
```

To setup a conda environment capable of building Larch, use:

`conda create -n larch` \
`conda activate larch` \
`conda install --channel "conda-forge" --update-deps --override-channels cmake make cxx-compiler openmpi openmpi-mpicc openmpi-mpicxx boost-cpp automake autoconf libtool yasm ucx zlib`

To setup a conda environment capable of building Larch including development tools, create `larch-dev` using the environment
file provided:

`conda env create -f environment.yml`

Building
--------

```shell
git submodule update --init --recursive
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
make -j16
```

Cmake build options:
  - add `-DCMAKE_CXX_CLANG_TIDY="clang-tidy"` to enable clang-tidy.
  - add `-DUSE_ASAN=yes` to enable asan and ubsan.

Running
-------

From the build directory:

```shell
ln -s ../data
./larch-test
```

Passing *nocatch* to the tests executable will allow exceptions to escape, which is useful for debugging. A gdb session can be started with `gdb --args build/larch-test nocatch`.

larch-test options:

- `nocatch` allows test exceptions to escape, which is useful for debugging. A gdb session can be started with `gdb --args build/larch-test nocatch`.
- `--list` produces a list of all available tests, along with an ID number.
- `--range` runs tests by ID with a string of comma-separated range or single ID arguments [e.g. 1-5,7,9,12-13].
- `-tag` excludes tests with a given tag.
- `+tag` includes tests with a given tag.
- For example, the `-tag "slow"` removes tests which require an long runtime to complete.

Third-party
-----------

* Lohmann, N. (2022). JSON for Modern C++ (Version 3.10.5) [Computer software]. https://github.com/nlohmann
* Eric Niebler. Range library for C++14/17/20. https://github.com/ericniebler/range-v3
