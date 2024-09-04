#!/bin/bash
#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_kmers.log
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH -n 12

kmer_start=$(date +%s)

print_elapsed_time() {
    local start_time=$1
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    local elapsed_minutes=$(elapsed/60)
    echo -e "Elapsed time: $elapsed minutes \n"
}

module load jellyfish-2.2.7-gcc-4.8.5-ocxxede

mkdir -p kmer_data

if [ -n "${ONT_READS}" ]; then
        start_time_10=$(date +%s)
        mkdir ./kmer_data/ont_kmers

        max_jobs=12
        job_count=0

        while read -r line; do
                jellyfish count -m 15 -s 10 -C ./ont_fastqs/"${line}".fastq -o ./kmer_data/ont_kmers/"${line}".jf &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < contig_names.txt

        wait

        echo "kmer counting for ont data done" >> progress.log
        print_elapsed_time $start_time_10 >> progress.log


        start_time_11=$(date +%s)
        while read -r line; do
                jellyfish histo ./kmer_data/ont_kmers/"${line}".jf > ./kmer_data/ont_kmers/"${line}"_hist.txt &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < contig_names.txt

        wait

        echo "kmer histograms for ont data done" >> progress.log

        rm ./kmer_data/ont_kmers/*.jf

        print_elapsed_time $start_time_11 >> progress.log
fi


if [ -n "${HIFI_READS}" ]; then
        start_time_12=$(date +%s)
        mkdir ./kmer_data/hifi_kmers

        max_jobs=12
        job_count=0

        while read -r line; do
                jellyfish count -m 15 -s 10 -C ./hifi_fastqs/"${line}".fastq -o ./kmer_data/hifi_kmers/"${line}".jf &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < contig_names.txt

        wait

        echo "kmer counting for hifi data done" >> progress.log
        print_elapsed_time $start_time_12 >> progress.log


        start_time_13=$(date +%s)
        while read -r line; do
                jellyfish histo ./kmer_data/hifi_kmers/"${line}".jf > ./kmer_data/hifi_kmers/"${line}"_hist.txt &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < contig_names.txt

        wait

        echo "kmer histograms for hifi data done" >> progress.log

        rm ./kmer_data/hifi_kmers/*.jf

        print_elapsed_time $start_time_13 >> progress.log
fi


if [ -n "${RNA_FORWARD}" ]; then
        start_time_14=$(date +%s)
        mkdir ./kmer_data/rna_kmers

        max_jobs=12
        job_count=0

        while read -r line; do
                jellyfish count -m 15 -s 10 -C ./rna_fastqs/"${line}".fastq -o ./kmer_data/rna_kmers/"${line}".jf &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < contig_names.txt

        wait

        echo "kmer counting for rna data done" >> progress.log
        print_elapsed_time $start_time_14 >> progress.log


        start_time_15=$(date +%s)
        while read -r line; do
                jellyfish histo ./kmer_data/rna_kmers/"${line}".jf > ./kmer_data/rna_kmers/"${line}"_hist.txt &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < contig_names.txt

        wait

        echo "kmer histograms for rna data done" >> progress.log

        rm ./kmer_data/rna_kmers/*.jf

        print_elapsed_time $start_time_15 >> progress.log
fi

if [ -n "${DNA_FORWARD}" ]; then
        start_time_20=$(date +%s)
        mkdir ./kmer_data/DNA_Illumina_kmers

        max_jobs=12
        job_count=0

        while read -r line; do
                jellyfish count -m 15 -s 10 -C ./DNA_Illumina_fastqs/"${line}".fastq -o ./kmer_data/DNA_Illumina_kmers/"${line}".jf &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < contig_names.txt

        wait

        echo "kmer counting for Illumina DNA data done" >> progress.log
        print_elapsed_time $start_time_20 >> progress.log


        start_time_21=$(date +%s)
        while read -r line; do
                jellyfish histo ./kmer_data/DNA_Illumina_kmers/"${line}".jf > ./kmer_data/DNA_Illumina_kmers/"${line}"_hist.txt &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < contig_names.txt

        wait

        echo "kmer histograms for Illumina DNA data done" >> progress.log

        rm ./kmer_data/DNA_Illumina_kmers/*.jf

        print_elapsed_time $start_time_20 >> progress.log
fi



kmer_end=$(date +%s)
kmer_runtime=$((kmer_end - kmer_start))
kmer_runtime_minutes=$((kmer_runtime / 60))
echo "kmer graphs completed in $kmer_runtime_minutes minutes" >> progress.log
