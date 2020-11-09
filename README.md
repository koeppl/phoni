# PHONI 
## - Practical Heuristic ON Incremental matching statistics computation

This framework supports the currently memory-friendliest way to compute the matching statistics of a pattern on highly-repetitive texts,
given that the input text is precomputed with the [https://github.com/maxrossi91/moni](MONI) index.

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
git clone --branch phoni https://github.com/maxrossi91/moni
```

### Compile

```console
mkdir build
cd build; cmake ..
make
```

In our experiments we used the file
[http://dolomit.cs.tu-dortmund.de/chr19.1000.fa.xz](chr19.1000.fa.xz)
as our text dataset,
and
[http://dolomit.cs.tu-dortmund.de/chr19.10.fa.xz](chr19.10.fa.xz) 
as our pattern dataset.
We have a shell script `benchmark.sh` for an automatic benchmark.
For this to work, some variables in it has to be set, as this project does not ship with the other matching statistic algorithms
such as [https://github.com/odenas/indexed_ms](msfast) or [https://github.com/apachecom/rrepair](rrepair), 
meaning it is necessary to download and compile those projects individually, and the set the corresponding variables in `benchmark.sh` manually.
Finally, the output of `benchmark.sh` can be processed by [https://github.com/koeppl/sqlplot](sqlplots) to generate the plots shown in the paper.

To compute the naive PHONI variant evaluated in the paper, simple exchange `lceToRBounded` with `lceToR_NaiveBounded` in the file `include/ms/phoni.hpp`.

