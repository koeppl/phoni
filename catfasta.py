#!/usr/bin/env python3

import fastaparser
import sys



with open(sys.argv[1]) as fasta_file:
	parser = fastaparser.Reader(fasta_file, parse_method='quick')
	for seq in parser:
		print(seq.sequence, end='')
