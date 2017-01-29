# KPC_grabber
This will extract KPC (and user definied flanking regions) from fasta sequences (genome assemblies/PacBio/etc.) in .gz format (how they are downloaded from genbank/refseq. In order for the KPC sequence to be extracted, it must meet certian criteria:

	- The flanking regions cannot extend over the start or end of the contig/seqeunce.
	- The sequence must not contain any N's
	- The KPC gene must be greater than 98% similar to the canonical KPC-2 gene sequence, and the alignment must stretch over 880 bp of the total 882 bp KPC sequence.

If all criteria are met, the program extracts the KPC sequence (and flanking regions). If the KPC gene is on the negative strand, the program with extract the reverse complement so that all extracted sequences are in the same orientation.


###Example Usage (extract KPC from all .gz files in directory with 1000bp flanking regions
> python KPC_grabber.py -f *.gz -o output.fasta -l 1000
