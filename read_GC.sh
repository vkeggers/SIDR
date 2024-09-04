#!/bin/bash
#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_read_GC.log
#SBATCH --mail-user=vegge003@fiu.edu


#calculate average GC content of the reads spanning the contig

read_GC_start=$(date +%s)

########################################################################################################################################################

#pacbio hifi
if [ -n "${HIFI_READS}" ]; then

#titles for columns in the output file
    echo -e "contig\tread_GC_hifi" > ./stats/read_GC_hifi.txt


#loop through each line in contig_names.txt
    while read -r line; do
        #use awk to print the second of four lines (2nd line is sequence) in the subset fastq | substitute any G or C or g or c to an empty character and
        #count the number of substitutions. Save this count to a variable called gc_number
        gc_number=$(awk 'NR%4 == 2 {print}' hifi_fastqs/"${line}".fastq | awk '{count += gsub(/[gcGC]/, "");} END {print count}')

        #use awk to print the second of four lines (2nd line is sequence) in the subset fastq | remove newline characters | count all characters and
        #save to a variable called total bases
        total_bases=$(awk 'NR%4 == 2 {print}' hifi_fastqs/"${line}".fastq | tr -d '\n' | wc -c)

        #pass the variables to awk, divide gc_number by total_bases and multiply by 100
        gc_percent=$(awk -v total="$total_bases" -v value="$gc_number" 'BEGIN { printf((value/total)*100); }')

        #print and append the values for the current contig to the output file
        echo -e "$line\t$gc_percent" >> ./stats/read_GC_hifi.txt

    done < contig_names.txt

fi


###################################################################################################################################################

#Oxford Nanopore
if [ -n "${ONT_READS}" ]; then

    echo -e "contig\tread_GC_ont" > ./stats/read_GC_ont.txt

    while read -r line; do
        gc_number=$(awk 'NR%4 == 2 {print}' ont_fastqs/"${line}".fastq | awk '{count += gsub(/[gcGC]/, "");} END {print count}')
        total_bases=$(awk 'NR%4 == 2 {print}' ont_fastqs/"${line}".fastq | tr -d '\n' | wc -c)
        gc_percent=$(awk -v total="$total_bases" -v value="$gc_number" 'BEGIN { printf((value/total)*100); }')
        echo -e "$line\t$gc_percent" >> ./stats/read_GC_ont.txt
    done < contig_names.txt

fi


##################################################################################################################################################

#RNA
if [ -n "${RNA_FORWARD}" ]; then

    echo -e "contig\tread_GC_rna" > ./stats/read_GC_rna.txt

    while read -r line; do
        gc_number=$(awk 'NR%4 == 2 {print}' rna_fastqs/"${line}".fastq | awk '{count += gsub(/[gcGC]/, "");} END {print count}')
        total_bases=$(awk 'NR%4 == 2 {print}' rna_fastqs/"${line}".fastq | tr -d '\n' | wc -c)
        gc_percent=$(awk -v total="$total_bases" -v value="$gc_number" 'BEGIN { printf((value/total)*100); }')
        echo -e "$line\t$gc_percent" >> ./stats/read_GC_rna.txt
    done < contig_names.txt

fi

###################################################################################################################################################

#DNA_Illumina
if [ -n "${DNA_FORWARD}" ]; then

    echo -e "contig\tread_GC_DNA_Illumina" > ./stats/read_GC_DNA_Illumina.txt

    while read -r line; do
        gc_number=$(awk 'NR%4 == 2 {print}' DNA_Illumina_fastqs/"${line}".fastq | awk '{count += gsub(/[gcGC]/, "");} END {print count}')
        total_bases=$(awk 'NR%4 == 2 {print}' DNA_Illumina_fastqs/"${line}".fastq | tr -d '\n' | wc -c)
        gc_percent=$(awk -v total="$total_bases" -v value="$gc_number" 'BEGIN { printf((value/total)*100); }')
        echo -e "$line\t$gc_percent" >> ./stats/read_GC_DNA_Illumina.txt
    done < contig_names.txt

fi

###################################################################################################################################################


read_GC_end=$(date +%s)
read_GC_runtime=$((read_GC_end - read_GC_start))
read_GC_runtime_minutes=$((read_GC_runtime / 60))
echo "read_GC stats completed in $read_GC_runtime_minutes minutes" >> progress.log
