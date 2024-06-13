# Getting started

Requirements
------------

* GCC 7.5
* cmake 3.16

For Ubuntu 18.04 LTS the following commands installs the requirements:

```shell
sudo apt install --no-install-recommends git cmake make g++ mpi-default-dev libprotobuf-dev libboost-dev libboost-program-options-dev libboost-filesystem-dev libboost-iostreams-dev libboost-date-time-dev protobuf-compiler automake autoconf libtool nasm
```

To get a recent cmake, download from `https://cmake.org/download/`, for example:

```shell
wget https://github.com/Kitware/CMake/releases/download/v3.23.1/cmake-3.23.1-linux-x86_64.tar.gz
```

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

```shell
conda create -n larch
conda activate larch
conda install --channel "conda-forge" --update-deps --override-channels cmake make cxx-compiler openmpi openmpi-mpicc openmpi-mpicxx boost-cpp automake autoconf libtool yasm ucx zlib
```

To setup a conda environment capable of building Larch including development tools, create `larch-dev` using the environment
file provided:

```shell
conda env create -f environment.yml
```

Building
--------

There are 4 executables that are built automatically as part of the larch package and provide various methods for exploring tree space and manipulating DAGs/trees:
- `larch-test` is the suite of tests used to validate the various routines.
- `larch-usher` takes an input tree/DAG and explores tree space through SPR moves.
- `merge` utility is used to manipulate(e.g. combine, prune)DAGs/trees.
- `dag2dot` utility writes a provided protobuf file in dot format for easy viewing.

Note: If you run against memory limitations during the cmake step, you can regulate number of parallel threads with `export CMAKE_NUM_THREADS="8"` (reduce number as necessary).

To build all from `larch/` directory, run:

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

### larch-test

From the `larch/build/` directory:
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

### file formats

For all tools in this suite, a number of file formats are supported for loading and storing MATs and MADAGs. When passing filepaths as arguments, the file format can be explicitly specified with `--input-format/--output-format` options.  Alternatively, the program can infer the file format when filepath contains a recognized file extension.

File format options:
- `MADAG dagbin` Supported as input and output. `*.dagbin` is the recognized extension.
- `MADAG protobuf` Supported as input and output. `*.pb_dag` is the recognized extension, or using `*.pb` WITHOUT a `--MAT-refseq-file` option.
- `MAT protobuf` Supported as input only. `*.pb_tree` is the recognized extension, or using `*.pb` WITH a `--MAT-refseq-file` option.
- `MADAG json` Supported as input only. `*.json_dag` or `*.json` is the recognized extension.

### larch-usher

From the `larch/build/` directory:
```shell
./larch-usher -i ../data/testcase/tree_1.pb.gz -o output_dag.pb -c 10
```
This command runs 10 iterations of larch-usher on the provided tree, and writes the final result to the file `output_dag.pb`

larch-usher options:
- `-i,--input` [REQUIRED] Filepath to the input tree/DAG (accepted file formats are: MADAG protobuf, MAT protobuf, JSON, Dagbin).
- `-o,--output` [REQUIRED] Filepath to the output tree/DAG (accepted file formats are: MADAG protobuf, Dagbin).
- `-c,--count` [Default: 1] Number of larch-usher iterations to run.
- `-r,--MAT-refseq-file` [REQUIRED if provided input file is a MAT protobuf] Filepath to json reference sequence.
- `-v,--VCF-input-file` Filepath to VCF containing ambiguous sequence data.
- `-l,--logpath` [Default: `optimization_log`] Filepath to write summary log.
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
- `--input-format` [Default: format inferred by file extension] Specify the format of the input file. Options are: (`dagbin`, `pb`, `dag-pb`, `tree-pb`, `json`, `dag-json`)
- `--output-format` [Default: format inferred by file extension] Specify the format of the output file. Options are: (`dagbin`, `pb`, `dag-pb`)

### merge

From the `larch/build/` directory:
```shell
./merge -i ../data/testcase/tree_1.pb.gz -i ../data/testcase/tree_2.pb.gz -d -o merged_trees.pb
```
This executable takes a list of protobuf files and merges the resulting DAGs together into one.

merge options:
- `-i,--input` Filepath to the input Tree/DAG (accepted file formats are: MADAG protobuf, MAT protobuf, JSON, Dagbin).
- `-o,--output` [Default: `merged.dagbin`] Filepath to the output Tree/DAG (accepted file formats are: MADAG protobuf, Dagbin).
- `-r,--MAT-refseq-file` [REQUIRED if input protobufs are MAT protobuf format] Filepath to json reference sequence.
- `-t,--trim` Trim output (Default trimming method is trim to best parsimony).
- `--rf` Trim output to minimize RF distance to the provided DAG file (Ignored if `-t` flag is not provided).
- `-s,--sample` Write a sampled single tree from DAG to file, rather than the whole DAG.
- `--input-format` [Default: format inferred by file extension] Specify the format of the input file(s). Options are: (`dagbin`, `pb`, `dag-pb`, `tree-pb`, `json`, `dag-json`)
- `--output-format` [Default: format inferred by file extension] Specify the format of the output file. Options are: (`dagbin`, `pb`, `dag-pb`)
- `--rf-format` [Default: format inferred by file extension] Specify the format of the RF file. Options are: (`dagbin`, `pb`, `dag-pb`, `tree-pb`, `json`, `dag-json`)

### dag2dot

From the `larch/build/` directory:
```shell
./dag2dot -i ../data/testcase/full_dag.pb
```
This command writes the provided DAG in dot format to stdout.

dag2dot options:
- `-i,--input` Filepath to the input Tree/DAG (accepted file formats are: MADAG protobuf, MAT protobuf, JSON, Dagbin).
- `-o,--output` [Default: DOT written to stdout] Filepath to the output DOT file.
- `--input-format` [Default: format inferred by file extension] Specify the format of the input file. Options are: (`dagbin`, `pb`, `dag-pb`, `tree-pb`, `json`, `dag-json`)
- `--dag/--tree` [REQUIRED if file extension is *.pb] Specify whether input file is a DAG or a Tree.


Third-party
-----------

* Lohmann, N. (2022). JSON for Modern C++ (Version 3.10.5) [Computer software]. https://github.com/nlohmann
* Eric Niebler. Range library for C++14/17/20. https://github.com/ericniebler/range-v3
