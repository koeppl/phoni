#!/usr/bin/env python3

import fastaparser
import sys



with open(sys.argv[1]) as fasta_file:
	parser = fastaparser.Reader(fasta_file, parse_method='quick')
	for seq in parser:
                print(seq.header)
                s = seq.sequence
                #_as_string()
                nlen = int((len(s)*int(sys.argv[2]))/100)
                print(s[:nlen])
