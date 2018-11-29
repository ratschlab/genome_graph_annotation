# Genome Graph Annotation Schemes

Sparse binary relation representations for genome graphs annotations

This repository implements the following schemes for representing graph annotation:
* Column-major compressed
* Row-major flat
* Rainbowfish
* BinRel-WT
* BRWT
* Multi-BRWT
* ...

As an underlying graph structure, the following representations are implemented:
1. Hash-based de Bruijn graph
2. Complete de Bruijn graph, taking constant space

## Comparison

The figures below show the final size of two compressed binary relations
* Kingsford with 3.7 bln rows and 2,652 columns, density ~0.3%
* Refseq (family) with 1 bln rows and 3,173 columns, density ~3.8%

| Method  | Kingsford, Gb  | RefSeq, Gb |
| ------------- | ---: | ---: |
| Column  | 36.56  | 80.18 |
| Flat  | 41.21  | 121.60 |
| [BinRel-WT](https://github.com/dieram3/binrel_wt/) | 49.57  | N/A |
| BinRel-WT (sdsl)  | 31.40  | 46.84 |
| Rainbowfish  | 19.22  | N/A |
| BRWT  | 12.97  | 51.82 |
| Multi-BRWT  | **9.19**  | **42.75** |

### Prerequisites
- cmake 3.6.1
- GNU GCC with C++17 (gcc-8 or higher) or LLVM Clang (clang-7 or higher)
- HTSlib
- boost
- folly (optional)

All can be installed with [brew](https://brew.sh) or [linuxbrew](https://linuxbrew.sh)

#### For compiling with GNU GCC:
```
brew install gcc autoconf automake libtool cmake make htslib
brew install --build-from-source boost
(optional) brew install --build-from-source double-conversion gflags glog lz4 snappy zstd folly
brew install gcc@8
```
Then set the environment variables accordingly:
```
echo "\
# Use gcc-8 with cmake
export CC=\"\$(which gcc-8)\"
export CXX=\"\$(which g++-8)\"
" >> $( [[ "$OSTYPE" == "darwin"* ]] && echo ~/.bash_profile || echo ~/.bashrc )
```

#### For compiling with LLVM Clang:
```
brew install llvm libomp autoconf automake libtool cmake make htslib boost folly
```
Then set the environment variables accordingly:
```
echo "\
# OpenMP
export LDFLAGS=\"\$LDFLAGS -L$(brew --prefix libomp)/lib\"
export CPPFLAGS=\"\$CPPFLAGS -I$(brew --prefix libomp)/include\"
# Clang C++ flags
export LDFLAGS=\"\$LDFLAGS -L$(brew --prefix llvm)/lib -Wl,-rpath,$(brew --prefix llvm)/lib\"
export CPPFLAGS=\"\$CPPFLAGS -I$(brew --prefix llvm)/include\"
export CXXFLAGS=\"\$CXXFLAGS -stdlib=libc++\"
# Path to Clang
export PATH=\"$(brew --prefix llvm)/bin:\$PATH\"
# Use Clang with cmake
export CC=\"\$(which clang)\"
export CXX=\"\$(which clang++)\"
" >> $( [[ "$OSTYPE" == "darwin"* ]] && echo ~/.bash_profile || echo ~/.bashrc )
```


### Compile
1. `git clone --recursive https://github.com/ratschlab/genome_graph_annotation`
2. make sure all submodules are downloaded: `git submodule update --init --recursive`
3. install third-party libraries from `external-libraries/` following the corresponding istructions  
or simply run the following script
```bash
git submodule update --init --recursive

pushd external-libraries/sdsl-lite
./install.sh $(pwd)
popd

pushd external-libraries/libmaus2
cmake -DCMAKE_INSTALL_PREFIX:PATH=$(pwd)
make -j $(($(getconf _NPROCESSORS_ONLN) - 1))
make install
popd
```

4. go to the **build** directory `mkdir -p build && cd build`
5. compile by `cmake .. && make -j $(($(getconf _NPROCESSORS_ONLN) - 1))`
6. run unit tests `./unit_tests`

#### Typical issues
* Linking against dynamic libraries in Anaconda when compiling libmaus2
  * make sure that packages like Anaconda are not listed in the exported environment variables

### Build types: `cmake .. <arguments>` where arguments are:
- `-DCMAKE_BUILD_TYPE=[Debug|Release|Profile]` -- build modes (`Release` by default)
- `-DBUILD_STATIC=[ON|OFF]` -- link statically (`OFF` by default)
- `-DWITH_AVX=[ON|OFF]` -- compile with support for the avx instructions (`ON` by default)

## Typical workflow
1. Build de Bruijn graph from Fasta files, FastQ files, or [KMC k-mer counters](https://github.com/refresh-bio/KMC/):\
`./annograph build`
2. Annotate graph using the column compressed annotation:\
`./annograph annotate`
3. Transform the built annotation to a different annotation scheme:\
`./annograph transform_anno`
4. Merge annotations (optional):\
`./annograph merge_anno`
5. Query annotated graph\
`./annograph classify`

### Example
```
./annograph build -k 12 -o tiny_example ../tests/data/tiny.fa

./annograph annotate -i tiny_example --anno-filename -o tiny_example ../tests/data/tiny.fa

./annograph classify -i tiny_example -a tiny_example.column.annodbg ../tests/data/tiny.fa

./annograph stats -a tiny_example tiny_example
```

## Scripts
For real benchmarking scripts, see [scripts](./scripts).

## Experiments on simulated matrices
Compressed simulated binary relation matrices can be generated using the script `experiments/run_benchmarks.py`. Given a column
count `$N_COLUMNS`, the simulation mode `$MODE` three available simulation modes are
1. `norepl` uniformly random matrix of size 1,000,000 x `$N_COLUMNS`,
2. `uniform_rows` 200,000 rows of size `$N_COLUMNS` duplicated 5 times to form a 1,000,000 x `$N_COLUMNS` matrix, and
3. `uniform_columns` `$N_COLUMNS` / 5 columns of size 1,000,000 duplicated 5 times to form a 1,000,000 x `$N_COLUMNS` matrix.

The compressor `$METHOD` can be one of: `brwt`, `bin_rel_wt`, `bin_rel_wt_sdsl`, `column`, `rbfish`, `flat`.

An experiment can then be run with the command
```
run_benchmarks.py $METHOD $MODE $N_COLUMNS $N_THREADS
```
when `.` is passed in place of `$METHOD` and/or `$MODE`, all methods/modes are run.

For method `brwt`, additional parameters can be passed at the end of the command. These can be one of
1. `--arity <N>` generate BRWT of arity `N`, or
2. `--greedy 1 --relax <N>` greedy optimization of column arrangement before construction of a BRWT of maximum arity `N`

All resulting matrices are saved to the `simulate` folder in the directory where the script is run.

### Figures from manuscript
To reproduce the simulated matrix experiment results from the [manuscript](https://www.biorxiv.org/content/early/2018/11/12/468512.full.pdf), run the following commands
```
for N_COLUMNS in 500 1000 3000; do
    run_benchmarks.py . . $N_COLUMNS $N_THREADS 
    run_benchmarks.py brwt $N_COLUMNS $N_THREADS --greedy 1 --relax 7
done
```

To plot all data for the figures in the experiments, run the command
```
run_benchmarks.py plot $N_COLUMNS
```
An alternative set of methods can be passed as subsequent arguments if desired, for example
```
run_benchmarks.py plot $N_COLUMNS brwt_arity_2 brwt_greedy_relax_6 bin_rel_wt
```

### Reference
> Mikhail Karasikov, Harun Mustafa, Amir Joudaki, Sara Javadzadeh No, Gunnar Rätsch, and André Kahles. [_Sparse Binary Relation Representations for Genome Graph Annotation_](https://www.biorxiv.org/content/early/2018/11/12/468512.full.pdf). (2018). https://doi.org/10.1101/468512.
