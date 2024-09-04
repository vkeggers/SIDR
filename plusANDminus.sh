#!/bin/bash
#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_read_GC.log
#SBATCH --mail-user=vegge003@fiu.edu


plusANDminus_start=$(date +%s)


#calculate the number of reads on the plus(forward) strand of DNA and minus(reverse) strand

#load modules
module load samtools-1.15.1-gcc-8.2.0


###########################################################################################################################################

#pacbio hifi
if [ -n "${HIFI_READS}" ]; then

        echo -e "contig\thifi_plus_counts" > ./stats/hifi_plus_counts.txt

        max_jobs=12
        job_count=0

        while read -r line; do

                samtools view -F 16 ./samsANDbams/hifi_alignments/"${line}"_reads.bam | awk -v line="$line" '{count++; contig=$3} END {if (NR==0) {print line 0} else {print line count}}' >> ./stats/hifi_plus_counts.txt &

                job_count=$((job_count + 1))

                        if [ "$job_count" -ge "$max_jobs" ]; then
                                wait
                                job_count=0
                        fi

        done < contig_names.txt

        wait

        sed -i 's/ /\t/' ./stats/hifi_plus_counts.txt


        echo -e "contig\thifi_minus_counts" > ./stats/hifi_minus_counts.txt

        while read -r line; do
                samtools view -f 16 ./samsANDbams/hifi_alignments/"${line}"_reads.bam | awk -v line="$line" '{count++; contig=$3} END {if (NR==0) {print line 0} else {print line count}}' >> ./stats/hifi_minus_counts.txt &

                job_count=$((job_count + 1))

                        if [ "$job_count" -ge "$max_jobs" ]; then
                                wait
                                job_count=0
                        fi

        done < contig_names.txt

        wait

        sed -i 's/ /\t/' ./stats/hifi_minus_counts.txt

fi


###########################################################################################################################################

#oxford nanopore
if [ -n "${ONT_READS}" ]; then

        echo -e "contig\tont_plus_counts" > ./stats/ont_plus_counts.txt

        max_jobs=12
        job_count=0

        while read -r line; do

                samtools view -F 16 ./samsANDbams/ont_alignments/"${line}"_reads.bam | awk -v line="$line" '{count++; contig=$3} END {if (NR==0) {print line 0} else {print line count}}' >> ./stats/ont_plus_counts.txt &

                job_count=$((job_count + 1))

                        if [ "$job_count" -ge "$max_jobs" ]; then
                                wait
                                job_count=0
                        fi

        done < contig_names.txt

        wait

        sed -i 's/ /\t/' ./stats/ont_plus_counts.txt


        echo -e "contig\tont_minus_counts" > ./stats/ont_minus_counts.txt

        while read -r line; do
                samtools view -f 16 ./samsANDbams/ont_alignments/"${line}"_reads.bam | awk -v line="$line" '{count++; contig=$3} END {if (NR==0) {print line 0} else {print line count}}' >> ./stats/ont_minus_counts.txt &

                job_count=$((job_count + 1))

                        if [ "$job_count" -ge "$max_jobs" ]; then
                                wait
                                job_count=0
                        fi

        done < contig_names.txt

        wait

        sed -i 's/ /\t/' ./stats/ont_minus_counts.txt

fi

###########################################################################################################################################

#RNA
if [ -n "${RNA_FORWARD}" ]; then

        echo -e "contig\trna_plus_counts" > ./stats/rna_plus_counts.txt

        max_jobs=12
        job_count=0

        while read -r line; do

                samtools view -F 16 ./samsANDbams/rna_alignments/"${line}"_reads.bam | awk -v line="$line" '{count++; contig=$3} END {if (NR==0) {print line 0} else {print line count}}' >> ./stats/rna_plus_counts.txt &

                job_count=$((job_count + 1))

                        if [ "$job_count" -ge "$max_jobs" ]; then
                                wait
                                job_count=0
                        fi

        done < contig_names.txt

        wait

        sed -i 's/ /\t/' ./stats/rna_plus_counts.txt


        echo -e "contig\trna_minus_counts" > ./stats/rna_minus_counts.txt

        while read -r line; do
                samtools view -f 16 ./samsANDbams/rna_alignments/"${line}"_reads.bam | awk -v line="$line" '{count++; contig=$3} END {if (NR==0) {print line 0} else {print line count}}' >> ./stats/rna_minus_counts.txt &

                job_count=$((job_count + 1))

                        if [ "$job_count" -ge "$max_jobs" ]; then
                                wait
                                job_count=0
                        fi

        done < contig_names.txt

        wait

        sed -i 's/ /\t/' ./stats/rna_minus_counts.txt

fi

############################################################################################################################################

#DNA Illumina
if [ -n "${DNA_FORWARD}" ]; then

        echo -e "contig\tDNA_illumina_plus_counts" > ./stats/DNA_Illumina_plus_counts.txt

        max_jobs=12
        job_count=0

        while read -r line; do

                samtools view -F 16 ./samsANDbams/DNA_Illumina_alignments/"${line}"_reads.bam | awk -v line="$line" '{count++; contig=$3} END {if (NR==0) {print line 0} else {print line count}}' >> ./stats/DNA_Illumina_plus_counts.txt &

                job_count=$((job_count + 1))

                        if [ "$job_count" -ge "$max_jobs" ]; then
                                wait
                                job_count=0
                        fi

        done < contig_names.txt

        wait

        sed -i 's/ /\t/' ./stats/DNA_Illumina_plus_counts.txt


        echo -e "contig\tDNA_Illumina_minus_counts" > ./stats/DNA_Illumina_minus_counts.txt

        while read -r line; do
                samtools view -f 16 ./samsANDbams/DNA_Illumina_alignments/"${line}"_reads.bam | awk -v line="$line" '{count++; contig=$3} END {if (NR==0) {print line 0} else {print line count}}' >> ./stats/DNA_Illumina_minus_counts.txt &

                job_count=$((job_count + 1))

                        if [ "$job_count" -ge "$max_jobs" ]; then
                                wait
                                job_count=0
                        fi

        done < contig_names.txt

        wait

        sed -i 's/ /\t/' ./stats/DNA_Illumina_minus_counts.txt

fi

############################################################################################################################################


plusANDminus_end=$(date +%s)
plusANDminus_runtime=$((plusANDminus_end - plusANDminus_start))
plusANDminus_runtime_minutes=$((plusANDminus_runtime / 60))
echo "plus and minus counts completed in $plusANDminus_runtime_minutes minutes" >> progress.log
