#!/bin/python

import subprocess, argparse, gzip, os
from Bio import SeqIO

#initiate parser
parser = argparse.ArgumentParser(description="This program will extract the KPC gene and flanking regions from fasta files")


#optional arguments
parser.add_argument("-f", "--files", help="file(s) to be analyzed", nargs="*", type=str)
parser.add_argument("-o", "--filename", help="name of output file to be created in current directory", default="results.fasta", type=str)
parser.add_argument("-l", "--length", help="length of flanking regions", default=0, type=int)

args = parser.parse_args()

#create file with KPC gene in it, create blast db with it
KPC_fasta = open("tmp.txt", "w")
KPC_fasta.write(">KPC2_882\nATGTCACTGTATCGCCGTCTAGTTCTGCTGTCTTGTCTCTCATGGCCGCTGGCTGGCTTTTCTGCCACCGCGCTGACCAACCTCGTCGCGGAACCATTCGCTAAACTCGAACAGGACTTTGGCGGCTCCATCGGTGTGTACGCGATGGATACCGGCTCAGGCGCAACTGTAAGTTACCGCGCTGAGGAGCGCTTCCCACTGTGCAGCTCATTCAAGGGCTTTCTTGCTGCCGCTGTGCTGGCTCGCAGCCAGCAGCAGGCCGGCTTGCTGGACACACCCATCCGTTACGGCAAAAATGCGCTGGTTCCGTGGTCACCCATCTCGGAAAAATATCTGACAACAGGCATGACGGTGGCGGAGCTGTCCGCGGCCGCCGTGCAATACAGTGATAACGCCGCCGCCAATTTGTTGCTGAAGGAGTTGGGCGGCCCGGCCGGGCTGACGGCCTTCATGCGCTCTATCGGCGATACCACGTTCCGTCTGGACCGCTGGGAGCTGGAGCTGAACTCCGCCATCCCAGGCGATGCGCGCGATACCTCATCGCCGCGCGCCGTGACGGAAAGCTTACAAAAACTGACACTGGGCTCTGCACTGGCTGCGCCGCAGCGGCAGCAGTTTGTTGATTGGCTAAAGGGAAACACGACCGGCAACCACCGCATCCGCGCGGCGGTGCCGGCAGACTGGGCAGTCGGAGACAAAACCGGAACCTGCGGAGTGTATGGCACGGCAAATGACTATGCCGTCGTCTGGCCCACTGGGCGCGCACCTATTGTGTTGGCCGTCTACACCCGGGCGCCTAACAAGGATGACAAGCACAGCGAGGCCGTCATCGCCGCTGCGGCTAGACTCGCGCTCGAGGGATTGGGCGTCAACGGGCAGTAA")
KPC_fasta.close()
subprocess.check_output("makeblastdb -in tmp.txt -dbtype nucl", shell=True)

for genome_file in args.files:

	run_blastn = subprocess.check_output("gzip -dc " + genome_file + " | blastn -db tmp.txt -outfmt 6", shell=True)
blast_results = run_blastn.rstrip().split('\n')

	#if no blast results, delete file
	if blast_results[0] == "":
		os.remove(genome_file)
		continue

	#Open file (append only)
	with open(args.filename, "a") as results_file:

		for result in blast_results:
			result = result.split('\t')

			#check here that % identity (98% or higher) and length (880bp or higher) of alignment meets standards
			if result[2] <= 98 or result[3] <= 880:
				continue

			#initialize blank variable to store results
			all_hits = ""
         
		        #if passes, extract full sequnce record using BioPython
			unzipped_genome_file = gzip.open(genome_file, 'rb')
			full_fasta = SeqIO.parse(unzipped_genome_file, "fasta")

			for seq in full_fasta:
			        if seq.id == result[0]:

				        #Check that KPC with flanking regions will not overlap with ends of contig, if so skip
					query_start = int(result[6]) - 1 #comply with python numbering
					query_end = int(result[7])
					flanking_region_size = args.length

					if query_start < flanking_region_size or query_end < flanking_region_size:
						continue
					elif query_start > (len(seq.seq) - flanking_region_size) or query_end > (len(seq.seq) - flanking_region_size):
						continue
            
					#Check that sequence doesn't contain and N's
					if "N" in seq.seq[query_start-flanking_region_size:query_end+flanking_region_size]:
						continue
					
					#If on - strand, take reverse complement
					finished_seq = seq.seq[query_start-flanking_region_size:query_end+flanking_region_size]
					if result[8] > result[9]:
						finished_seq = finished_seq.reverse_complement()

					#If all criteria passes, append the result to all_hits
					all_hits = all_hits + ">" + str(seq.id) + "\n" + finished_seq + "\n"

			results_file.write(str(all_hits))

#remove tmp files
os.remove("tmp.txt")
os.remove("tmp.txt.nhr")
os.remove("tmp.txt.nin")
os.remove("tmp.txt.nsq")
