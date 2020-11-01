#!/usr/bin/env python3
from Bio import SeqIO
import sys

i=0
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
	with open(sys.argv[2] + "/" + str(i), 'w') as file:
		file.write(str(seq_record.seq))
	i+=1

