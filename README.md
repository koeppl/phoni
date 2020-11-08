# moni
A Read Aligner with Multi-Genome References


pip3 install biopython
pip3 install fastaparser
pip3 install psutil
ln -sv build/no_thresholds .
wget 'http://dolomit.cs.tu-dortmund.de/chr19.10.fa.xz'; wget 'http://dolomit.cs.tu-dortmund.de/chr19.1000.fa.xz'


for phoni_0, you have to exchange lceToRBounded with lceToR_NaiveBounded in include/ms/phoni.hpp
