#!/bin/bash
#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_depth.log
#SBATCH --mail-user=vegge003@fiu.edu

depth_start=$(date +%s)

module load samtools-1.15.1-gcc-8.2.0

if [ -n "${HIFI_READS}" ]; then

        echo -e "contig\tAvg_depth_hifi" > ./stats/Avg_depth_hifi.txt

        while read -r line; do

                Avg_depth_hifi=$(samtools depth ./samsANDbams/hifi_alignments/"${line}"_reads.bam | \
                awk '{sum+=$3; count++} END {if (count > 0) print sum/count; else print "0"}')
                echo -e "${line}\t${Avg_depth_hifi}" >> ./stats/Avg_depth_hifi.txt

        done < contig_names.txt

        wait

        echo "Average depth of hifi reads finished" >> progress.log

fi


if [ -n "${ONT_READS}" ]; then

        echo -e "contig\tAvg_depth_ont" > ./stats/Avg_depth_ont.txt

        while read -r line; do

                Avg_depth_ont=$(samtools depth ./samsANDbams/ont_alignments/"${line}"_reads.bam | \
                awk '{sum+=$3; count++} END {if (count > 0) print sum/count; else print "0"}')
                echo -e "${line}\t${Avg_depth_ont}" >> ./stats/Avg_depth_ont.txt

        done < contig_names.txt

        wait

        echo "Average depth of ont reads finished" >> progress.log

fi


if [ -n "${RNA_FORWARD}" ]; then

        echo -e "contig\tAvg_depth_rna" > ./stats/Avg_depth_rna.txt

        while read -r line; do

                Avg_depth_rna=$(samtools depth ./samsANDbams/rna_alignments/"${line}"_reads.bam | \
                awk '{sum+=$3; count++} END {if (count > 0) print sum/count; else print "0"}')
                echo -e "${line}\t${Avg_depth_rna}" >> ./stats/Avg_depth_rna.txt

        done < contig_names.txt

        wait

        echo "Average depth of rna reads finished" >> progress.log

fi


if [ -n "${DNA_FORWARD}" ]; then

        echo -e "contig\tAvg_depth_DNA_Illumina" > ./stats/Avg_depth_DNA_Illumina.txt

        while read -r line; do

                Avg_depth_DNA_Illumina=$(samtools depth ./samsANDbams/DNA_Illumina_alignments/"${line}"_reads.bam | \
                awk '{sum+=$3; count++} END {if (count > 0) print sum/count; else print "0"}')
                echo -e "${line}\t${Avg_depth_DNA_Illumina}" >> ./stats/Avg_depth_DNA_Illumina.txt

        done < contig_names.txt

        wait

        echo "Average depth of Illumina DNA reads finished" >> progress.log

fi

depth_end=$(date +%s)
depth_runtime=$((depth_end - depth_start))
depth_runtime_minutes=$((depth_runtime / 60))
echo "Depth stats completed in $depth_runtime_minutes minutes" >> progress.log
