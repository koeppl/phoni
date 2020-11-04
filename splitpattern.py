#!/usr/bin/env python3
from Bio import SeqIO
import sys

descriptions=[]
i=0
for seq_record in SeqIO.parse(sys.argv[1], "fasta"):
	descriptions.append(seq_record.description)
	with open(sys.argv[2] + "/" + str(i), 'w') as file:
		file.write(str(seq_record.seq))
	i+=1

with open(sys.argv[2] + "/desc.txt" , 'w') as f:
	for desc in descriptions:
		print(desc, file=f)
