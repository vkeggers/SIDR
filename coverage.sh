#!/bin/bash
#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_coverage.log
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH -n 12

#various forms of coverage


coverage_start=$(date +%s)


#########################################################################################################################################

#pacbio hifi

if [ -n "${HIFI_READS}" ]; then

#titles for columns in the output file
echo -e "contig\thifi_average_coverage\thifi_coverage_percent" > ./stats/hifi_coverageStats.txt

#samtools depth will give you the number of reads that cover each nucleotide position in the genome
samtools depth ./samsANDbams/hifi_alignments/HIFIaln_sorted.bam > hifi_coverage.txt

#If column 3 of the samtools depth output (coverage at that nucleotide) is greater than zero, add 1 to covered_bases and add the value to the column to sum_
#coverage. If column 3 is not greater than zero, just add 1 to total_bases, which is basically just the length of the contig
#The END statement is a for loop that calculates and prints the average coverage over the whole contig and the coverage percent (telling you if any bases we
#re not covered by the reads) for each contig. These values are appended to the output file (PBcoverageStats.txt) which we initialized earlier
awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum_coverage[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum_coverage) {
        average_coverage = sum_coverage[contig] / total_bases[contig];
        coverage_percent = ((covered_bases[contig] / total_bases[contig]) *100);
        print contig, average_coverage, coverage_percent;
    }
}' hifi_coverage.txt >> ./stats/hifi_coverageStats.txt

#replaces spacecs with tabs, will be important when concatenating the tables together at the end
sed -i 's/ /\t/g' ./stats/hifi_coverageStats.txt

fi

###################################################################################################################################

#oxford nanopore

if [ -n "${ONT_READS}" ]; then

echo -e "contig\tont_average_coverage\tont_coverage_percent" > ./stats/ont_coverageStats.txt
samtools depth ./samsANDbams/ont_alignments/ONTaln_sorted.bam > ont_coverage.txt

awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum_coverage[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum_coverage) {
        average_coverage = sum_coverage[contig] / total_bases[contig];
        coverage_percent = ((covered_bases[contig] / total_bases[contig]) *100);
        print contig, average_coverage, coverage_percent;
    }
}' ont_coverage.txt >> ./stats/ont_coverageStats.txt
sed -i 's/ /\t/g' ./stats/ont_coverageStats.txt

fi

####################################################################################################################################

#RNA

if [ -n "${RNA_FORWARD}" ]; then

echo -e "contig\trna_average_coverage\trna_coverage_percent" > ./stats/rna_coverageStats.txt
samtools depth ./samsANDbams/rna_alignments/RNAaln_sorted.bam > rna_coverage.txt

awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum_coverage[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum_coverage) {
        average_coverage = sum_coverage[contig] / total_bases[contig];
        coverage_percent = ((covered_bases[contig] / total_bases[contig]) *100);
        print contig, average_coverage, coverage_percent;
    }
}' rna_coverage.txt >> ./stats/rna_coverageStats.txt
sed -i 's/ /\t/g' ./stats/rna_coverageStats.txt

fi

####################################################################################################################################

#DNA_Illumina

if [ -n "${DNA_FORWARD}" ]; then

echo -e "contig\tDNA_Illumina_average_coverage\tDNA_Illumina_coverage_percent" > ./stats/DNA_Illumina_coverageStats.txt
samtools depth ./samsANDbams/DNA_Illumina_alignments/DNA_Illumina_aln_sorted.bam > DNA_Illumina_coverage.txt

awk '{
    if ($3 > 0) {
        covered_bases[$1]++;
        sum_coverage[$1] += $3;
    }
    total_bases[$1]++;
}
END {
    for (contig in sum_coverage) {
        average_coverage = sum_coverage[contig] / total_bases[contig];
        coverage_percent = ((covered_bases[contig] / total_bases[contig]) *100);
        print contig, average_coverage, coverage_percent;
    }
}' DNA_Illumina_coverage.txt >> ./stats/DNA_Illumina_coverageStats.txt
sed -i 's/ /\t/g' ./stats/DNA_Illumina_coverageStats.txt

fi

####################################################################################################################################

coverage_end=$(date +%s)
coverage_runtime=$((coverage_end - coverage_start))
coverage_runtime_minutes=$((coverage_runtime / 60))
echo "coverage stats completed in $coverage_runtime_minutes minutes" >> progress.log
