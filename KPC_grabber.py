#!/bin/python

import subprocess
from Bio import SeqIO

run_blastn = subprocess.check_output("blastn -db /home/johnsonrc2/work/KPC_grabber/KPC.fasta -query /home/johnsonrc2/work/KPC_grabber/test_genomes/$
blast_results = run_blastn.rstrip().split('\n')

for result in blast_results:
        result = result.split('\t')
        full_fasta = SeqIO.parse("/home/johnsonrc2/work/KPC_grabber/test_genomes/multi_KPC/Kpne_PLA68_122414_abfpw.masurca.fasta", "fasta")
        for seq in full_fasta:
                if seq.id == result[0]:
                        print(">" + str(seq.id))     
                        print(seq.seq[int(result[6]):int(result[7])])
