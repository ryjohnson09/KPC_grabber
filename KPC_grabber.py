#!/bin/python

import subprocess
from Bio import SeqIO

run_blastn = subprocess.check_output("blastn -db /home/johnsonrc2/work/KPC_grabber/KPC.fasta -query /home/johnsonrc2/work/KPC_grabber/test_genomes/$
blast_results = run_blastn.rstrip().split('\n')

#Open file (append only)               
with open("KPC_results.fasta", "a") as results_file:

        for result in blast_results:
                result = result.split('\t')
                #check here that % identity (98% or higher) and length (880bp or higher) of alignment meets standards
                if float(result[2]) <= 98 or float(result[3]) <= 880:
                        continue

                #initialize blank variable to store results
                all_hits = ""

                #if passes, extract full sequnce record using BioPython
                full_fasta = SeqIO.parse("/home/johnsonrc2/work/KPC_grabber/test_genomes/multi_KPC/Kpne_PLA68_122414_abfpw.masurca.fasta", "fasta")
                for seq in full_fasta:
                        if seq.id == result[0]:

                                #Check that KPC with flanking regions will not overlap with ends of contig, if so skip
                                query_start = int(result[6])
                                query_end = int(result[7])
                                flanking_region_size = 1000
                                if query_start < flanking_region_size or query_end < flanking_region_size:
                                        continue
                                elif query_start > (len(seq.seq) - flanking_region_size) or query_end > (len(seq.seq) - flanking_region_size):
                                        continue

                                #If all criteria passes, append the result to all_hits
                                all_hits = all_hits + ">" + str(seq.id) + "\n" + seq.seq[query_start:query_end] + "\n"

                results_file.write(str(all_hits))
