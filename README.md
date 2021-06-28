# PHONI 
(Practical Heuristic ON Incremental matching statistics computation)

This framework supports the currently memory-friendliest way to compute the matching statistics of a pattern on highly-repetitive texts,
given that the input text is precomputed with the [MONI](https://github.com/maxrossi91/moni) index 
(more precisely, we need all ingredients of MONI except the thresholds).

We require the pattern and the text to be available in form of sequences stored in the `.fa` (FASTA) format.
To use our solution, you need to have recent `cmake`, `g++`, `zsh`, and `python 3` installed.

## Preparations

We need the following python 3 packages for extracting and concatenating `.fa` files:
```console
	pip3 install biopython
	pip3 install fastaparser
	pip3 install psutil
```


```console
git clone --branch phoni https://github.com/koeppl/phoni
```

### Compile

```console
mkdir build
cd build; cmake ..
make
```

### Building the index

To *build* the index we use the command `phoni build` from the build directory.

``` console
phoni build \
-r <filename of the reference> \
-t <number of threads> \
-g <grammar format> \
-f <input file is a fasta file> \
```
For example, to build the phoni index for the file `yeast.fasta` using 4 `threads` and the `plain` grammar we run from the `build` folder:
``` conole
python3 phoni build -r ../data/yeast.fasta -f -t 4 -g plain
```

This command will produce `yeast.fasta.phoni` and `yeast.fasta.plain.slp` in the `data` folder, which represent the `phoni` index.

### Querying the index

To *query* the index we use the command `phoni query` from the build directory.

``` console
phoni ms \
-r <filename of the reference> \
-p <number of threads> \
-g <grammar format> \
```
For example, to query the phoni index for the file `yeast.fasta` using the `plain` grammar with the pattern `samples.fastq` we run from the `build` folder:
``` conole
python3 phoni build -r ../data/yeast.fasta -p ../data/samples.fa -g plain
```

This command will produce `samples.fa.positions` and `samples.fa.lengths` in the `data` folder, which represent the matching staistics *positions* and *lengths* of `samples.fa` against `yeast.fasta`, respectively.

### Compatibility queries

To perform the queries using the old *query* command we first have to split the `samples.fa` file using the python tool `splitpattern.py`:

```console
mkdir data/samples.fa.dir
python3 splitpattern.py data/samples.fa data/samples.fa.dir
```

Then we can query `samples.fa` against the `yeast.fasta` index with the plain grammar using the following command:

```console
./build/test/src/phoni_compatibility data/yeast.fasta -p data/samples.fa -g plain
```

<!-- ### Run

After building, the `build` directory should contain several python scripts for building the index data structures:

 - `moni` for building all MONI data structures 
 - `no_thresholds` : `moni` except computing the thresholds
 - `thresholds`: compute only the thresholds

The syntax of all these scripts is the same as described in [MONI](https://github.com/maxrossi91/moni).
All these scripts create auxiliary files of a given text `text.fa` whose filenames consists of an additional file extension of `text.fa`.

Finally, to run `phoni`, we first build the RLBWT and its auxiliary data structures with `./test/src/build_phoni -f text.fa`.
Then we can run `./test/src/phoni -f text.fa -p pattern.fa` to compute the matching statistics of `pattern.fa` within the text `text.fa`. -->

### python Tools

 - `prefixpattern.py`: takes the x% prefix of each pattern stored in a `.fa` file and outputs a new `.fa` file to stdout
 - `splitpattern.py`: splits a `.fa` file into individual sequences stored in a directory given as program parameter (we need this for `msfast` as it does not support reading `.fa` files)

### Benchmarks

We provide a script and benchmark files to evaluate PHONI in the setting as described in the paper

Christina Boucher, Travis Gagie, Tomohiro I, Dominik Köppl, Ben Langmead, Giovanni Manzini, Gonzalo Navarro, Alejandro Pacheco, Massimiliano Rossi: [PHONI: Streamed Matching Statistics with Multi-Genome References](https://doi.org/10.1109/DCC50243.2021.00027), In proceedings of 2021 Data Compression Conference, pp. 193-202, 2021. 

Christina Boucher, Travis Gagie, Tomohiro I, Dominik Köppl, Ben Langmead, Giovanni Manzini, Gonzalo Navarro, Alejandro Pacheco, Massimiliano Rossi: PHONI: Streamed Matching Statistics with Multi-Genome References, [arXiv:2011.05610](https://arxiv.org/abs/2011.05610), 11 Nov 2020.

In our experiments we used the file

 - [chr19.1000.fa.xz](http://dolomit.cs.tu-dortmund.de/chr19.1000.fa.xz) as our text dataset, and
 - [chr19.10.fa.xz](http://dolomit.cs.tu-dortmund.de/chr19.10.fa.xz) as our pattern dataset.

We have a shell script `benchmark.sh` for an automatic benchmark.
For this to work, some variables in it has to be set, as this project does not ship with the other matching statistic algorithms, namely

 - [MONI](https://github.com/maxrossi91/moni)
 - [msfast](https://github.com/odenas/indexed_ms), and
 - [rrepair](https://github.com/apachecom/rrepair).

meaning it is necessary to download and compile those projects individually, and the set the corresponding variables in `benchmark.sh` manually
(more precisely: in the switch-case statement for the hostname in the beginning).
Finally, the output of `benchmark.sh` can be processed by [sqlplots](https://github.com/koeppl/sqlplot) to generate the plots shown in the paper.

To compute the naive PHONI variant evaluated in the paper, simple exchange `lceToRBounded` with `lceToR_NaiveBounded` in the file `include/ms/phoni.hpp`.
