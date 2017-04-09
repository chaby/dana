This is code for a pipeline to take paired-end Illumina reads of amplicons and assemble them.

The code does the following:

INPUT: fastq files

1. Quality trim reads
2. Demultiplex the data and remove the primers
3. Choose sequence length for each marker (forward and reverse reads separately; 
	done by the investigator based on plots generaed using R scripts length_histograms.R and threshold_identification.R)
4. Trim sequences to desired size
5. Attach F and R reads (taking the reverse complement of R); if F& R reads do not overlap by at least 10 basepairs, reads are joined end-to-end
6. Remove singletons
7. Find unique sequences and count their number of reads

OUTPUT: fasta file of unique sequences, table of the counts of these sequences

Then, using R scripts,

1. For each individual, identify the most abundant read in the count table and compare it to the second most abundant read. 
	Select which individuals to keep. Done with counttables_DANAcode.R and count_table_DANA_diploid.R
2. extract the read(s) for each individual that we want to keep.
	The result is a fasta file where each individual is represented by 1 sequence (for mt loci) or 2 sequences (for diploid loci).
	This is done with extract_mito_DANAcode.R and extract_all_seqs_DIPLOID.R
