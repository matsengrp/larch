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

`singularity build larch-singularity.sif larch-singularity.def`
`singularity shell larch-singularity.sif --net`

To setup a conda environment capable of building Larch, use:

`conda create -n larch`
`conda activate larch`
`conda install --channel "conda-forge" --update-deps --override-channels cmake make cxx-compiler openmpi openmpi-mpicc openmpi-mpicxx boost-cpp automake autoconf libtool yasm ucx zlib`

To setup a conda environment capable of building Larch including development tools, create `larch-dev` using the environment
file provided:

`conda env create -f environment.yml`

Building
--------

There are 4 executables that are built automatically as part of the larch package and provide various methods for exploring tree space and manipulating DAGs/trees: 
- `larch-test` is the suite of tests used to validate the various routines.
- `larch-usher` takes an input tree/DAG and explores tree space through SPR moves.
- `merge` utility is used to manipulate(e.g. combine, prune)DAGs/trees.
- `dag2dot` utility writes a provided protobuf file in dot format for easy viewing.

Note: If you run against memory limitations during the cmake step, you can regulate number of parallel threads with `export CMAKE_NUM_THREADS="8"` (reduce number as necessary).

To build all from `larch/` directory, run:

`git submodule update --init --recursive`
`mkdir build`
`cd build`
`cmake -DCMAKE_BUILD_TYPE=Debug ..`
`make -j16`

Cmake build options:
  - add `-DCMAKE_CXX_CLANG_TIDY="clang-tidy"` to enable clang-tidy.
  - add `-DUSE_ASAN=yes` to enable asan and ubsan.

Running
-------

### larch-test

From the `larch/build/` directory:

`ln -s ../data`
`./larch-test`

Passing *nocatch* to the tests executable will allow exceptions to escape, which is useful for debugging. A gdb session can be started with `gdb --args build/larch-test nocatch`.

larch-test options:

- `nocatch` allows test exceptions to escape, which is useful for debugging. A gdb session can be started with `gdb --args build/larch-test nocatch`.
- `--list` produces a list of all available tests, along with an ID number.
- `--range` runs tests by ID with a string of comma-separated range or single ID arguments [e.g. 1-5,7,9,12-13].
- `-tag` excludes tests with a given tag.
- `+tag` includes tests with a given tag.
- For example, the `-tag "slow"` removes tests which require an long runtime to complete.

### larch-usher

From the `larch/build/` directory:
```shell
./larch-usher -i ../data/testcase/tree_1.pb.gz -o output_dag.pb -c 10
```
This command runs 10 iterations of larch-usher on the provided tree, and writes the final result to the file `output_dag.pb`

larch-usher options:
- `-i,--input` [REQUIRED] The name of the input tree/DAG (accepted file formats are: MADAG protobuf, MAT protobuf, JSON).
- `-o,--output` [REQUIRED] The file path to write the resulting DAG to.
- `-c,--count` [Default: 1] Number of larch-usher iterations to run.
- `-r,--MAT-refseq-file` [REQUIRED if provided input file is a MAT protobuf] Reference sequence file.
- `-v,--VCF-input-file` VCF file containing ambiguous sequence data.
- `-l,--logpath` [Default: `optimization_log`] Filepath to write log to.
- `-s,--switch-subtrees` [Default: never] Switch to optimizing subtrees after the specified number of iterations.
- `--min-subtree-clade-size` [Default: 100] The minimum number of leaves in a subtree sampled for optimization (ignored without option `-s`).
- `--max-subtree-clade-size` [Default: 1000] The maximum number of leaves in a subtree sampled for optimization (ignored without option `-s`).
- `--move-coeff-nodes` [Default: 1] New node coefficient for scoring moves. Set to 0 to apply only parsimony-optimal SPR moves.
- `--move-coeff-pscore` [Default: 1] Parsimony score coefficient for scoring moves. Set to 0 to apply only topologically novel SPR moves.
- `--sample-method` [Default: `parsimony`] Select method for sampling optimization tree from the DAG. Options are: (`parsimony`, `random`, `rf-minsum`, `rf-maxsum`).
- `--sample-uniformly` [Default: use natural distribution] Use a uniform distribution to sample trees for optimization.
- For example, if the sampling method is `parsimony` and `--sample-uniformly` is provided, then a uniform distribution on parsimony-optimal trees is sampled from.
- `--callback-option` [Default: `best-moves`] Specify which SPR moves are chosen and applied. Options are: (`all-moves`, `best-moves-fixed-tree`, `best-moves-treebased`, `best-moves`).
- `--trim` [Default: do not trim] Trim optimized dag to contain only parsimony-optimal trees before writing to protobuf.
- `--keep-fragment-uncollapsed` [Default: collapse] Do not collapse empty (non-mutation-bearing) edges in the optimization tree.
- `--quiet` [Default: write intermediate files] Do not write intermediate protobuf file at each iteration.

### merge

From the `larch/build/` directory:
```shell
./merge -i ../data/testcase/tree1.pb.gz -i ../data/testcase/tree2.pb.gz -d -o merged_trees.pb
```
This executable takes a list of protobuf files and merges the resulting DAGs together into one.

merge options:
- `-i,--input` Input protobuf files.
- `-o,--output` [Default: `merged.pb`] Save the output to filename.
- `-r,--refseq` [REQUIRED if input protobufs are MAT protobuf format] Read reference sequence from file.
- `-d,--dag` Input files are MADAG protobuf format\n";
- `-t,--trim` Trim output (default trimming method is trim to best parsimony).
- `--rf` Trim output to minimize RF distance to the provided protobuf(Ignored if `-t` flag is not provided).
- `-s,--sample` Write a sampled single tree from DAG to file, rather than the whole DAG.

### dag2dot

From the `larch/build/` directory:
```shell
./dag2dot -d ../data/testcase/full_dag.pb
```
This command writes the provided DAG in dot format to stdout.

dag2dot options:
- `-t,--tree-pb` Input MAT protobuf filename.
- `-d,--dag-pb` Input DAG protobuf filename.
- `-j,--dag-json` Input DAG json filename.


Third-party
-----------

* Lohmann, N. (2022). JSON for Modern C++ (Version 3.10.5) [Computer software]. https://github.com/nlohmann
* Eric Niebler. Range library for C++14/17/20. https://github.com/ericniebler/range-v3
