#!/usr/bin/env python3
from Bio import SeqIO
import sys

s=''
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
	# print(seq_record.id)
	# print(repr(seq_record.seq))
	# print(len(seq_record))
	s += seq_record.seq
print(s)

