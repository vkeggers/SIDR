# SIDR
### A genome decontamination software

The SIDR.sh file ties all the smaller scripts together. The smaller scripts contain parts of the program meant to calculate statistics for each contig.

The python file contains the machine learning program and is run after all genome stats have been calculated and concatenated into a table.

The test.fa is a short fasta file meant to ensure that statistics like GC content and sequence length are being calculated correctly. Similarly, the testAnswers.txt gives the true statistics for the sequences in test.fa.
