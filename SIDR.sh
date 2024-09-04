#!/bin/bash
#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_alignments.log
#SBATCH --mail-user=vegge003@fiu.edu
#SBATCH -n 12

export total_start=$(date +%s)
touch progress.log

export GENOME=/scratch/jfierst/tori/SIDR/testData/pooja/ONTtestData02_flye_assembly_output.fasta
export ONT_READS=/scratch/jfierst/tori/SIDR/testData/pooja/ONTtestData02.fastq
export RNA_FORWARD=/scratch/jfierst/tori/SIDR/testData/pooja/Illumina_RNA_testData_F.fastq
export RNA_REVERSE=/scratch/jfierst/tori/SIDR/testData/pooja/Illumina_RNA_testData_R.fastq
export DNA_FORWARD=
export DNA_REVERSE=
export HIFI_READS=
export NT=/home/data/jfierst/veggers/nt_db/nt
export SAMPLE=testData

module load samtools-1.15.1-gcc-8.2.0
module load blast-plus-2.7.1-gcc-8.2.0-fm5yf3k
module load minimap2-2.24
module load bwa-0.7.17-gcc-8.2.0-qgdird7
module load bedtools2-2.27.1-gcc-8.2.0-bxmhnwb

mkdir -p stats

print_elapsed_time() {
    local start_time=$1
    local end_time=$(date +%s)
    local elapsed=$(( end_time - start_time ))
    local elapsed_minutes=$(elapsed/60)
    echo -e "Elapsed time: $elapsed minutes \n"
}


samsANDbams_start=$(date +%s)

grep ">" $GENOME | sed 's/>//g' > contig_names.txt

mkdir contigs
awk '/^>/{filename=sprintf("%s/%s.fasta", "contigs", substr($0,2)); next} {print > filename}' $GENOME

mkdir samsANDbams

if [ -n "${HIFI_READS}" ]; then

        start_time_1=$(date +%s)

        mkdir hifi_fastqs
        cd samsANDbams
        mkdir hifi_alignments
        cd hifi_alignments

        minimap2 -ax map-hifi ${GENOME} ${HIFI_READS} -o HIFIaln.sam
        samtools view -Sb ./HIFIaln.sam -o ./HIFIaln.bam
        samtools sort -o ./HIFIaln_sorted.bam ./HIFIaln.bam
        samtools index ./HIFIaln_sorted.bam ./HIFIaln_sorted.bai
        bedtools bamtobed -i ./HIFIaln_sorted.bam > output.bed

        max_jobs=12
        job_count=0

        while read -r line; do
                grep "${line}" output.bed > "${line}"_output.bed &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt
        wait

        echo "Bed files for hifi reads done" >> ./../../progress.log
        print_elapsed_time $start_time_1 >> ./../../progress.log


        start_time_2=$(date +%s)
        while read -r line; do

                samtools view -L "${line}"_output.bed -h HIFIaln_sorted.bam | samtools view -b - > "${line}"_reads.bam &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait
        echo "Bam files for hifi reads done" >> ./../../progress.log
        print_elapsed_time $start_time_2 >> ./../../progress.log

        rm *.bed


        start_time_3=$(date +%s)
        while read -r line; do

                bedtools bamtofastq -i "${line}"_reads.bam -fq ./../../hifi_fastqs/"${line}".fastq &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Fastq files for hifi reads done" >> ./../../progress.log
        print_elapsed_time $start_time_3 >> ./../../progress.log

        cd ..
        cd ..
fi

wait

if [ -n "${ONT_READS}" ]; then

        start_time_4=$(date +%s)

        mkdir ont_fastqs
        cd samsANDbams
        mkdir ont_alignments
        cd ont_alignments

        minimap2 -ax map-ont ${GENOME} ${ONT_READS} -o ONTaln.sam
        samtools view -Sb ./ONTaln.sam -o ./ONTaln.bam
        samtools sort -o ./ONTaln_sorted.bam ./ONTaln.bam
        samtools index ./ONTaln_sorted.bam ./ONTaln_sorted.bai
        bedtools bamtobed -i ./ONTaln_sorted.bam > output.bed

        max_jobs=12
        job_count=0
         while read -r line; do
                grep "${line}" output.bed > "${line}"_output.bed &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Bed files for ont reads done" >> ./../../progress.log
        print_elapsed_time $start_time_4 >> ./../../progress.log

        start_time_5=$(date +%s)
        while read -r line; do
                samtools view -L "${line}"_output.bed -h ONTaln_sorted.bam | samtools view -b - > "${line}"_reads.bam &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Bam files for ont reads done" >> ./../../progress.log
        print_elapsed_time $start_time_5 >> ./../../progress.log

        rm *.bed


        start_time_6=$(date +%s)
        while read -r line; do
                bedtools bamtofastq -i "${line}"_reads.bam -fq ./../../ont_fastqs/"${line}".fastq &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Fastq files for ont reads done" >> ./../../progress.log
        print_elapsed_time $start_time_6 >> ./../../progress.log

        cd ..
        cd ..
fi

wait

if [ -n "${RNA_FORWARD}" ]; then

        start_time_7=$(date +%s)

        mkdir rna_fastqs
        cd samsANDbams
        mkdir rna_alignments
        cd rna_alignments

        bwa index ${GENOME}
        bwa mem -t 4 ${GENOME} ${RNA_FORWARD} ${RNA_REVERSE} > RNAaln.sam
        samtools view -Sb ./RNAaln.sam -o ./RNAaln.bam
        samtools sort -o ./RNAaln_sorted.bam ./RNAaln.bam
        samtools index ./RNAaln_sorted.bam ./RNAaln_sorted.bai
        bedtools bamtobed -i ./RNAaln_sorted.bam > output.bed

        max_jobs=12
        job_count=0

        while read -r line; do
                grep "${line}" output.bed > "${line}"_output.bed &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Bed files for rna reads done" >> ./../../progress.log
        print_elapsed_time $start_time_7 >> ./../../progress.log

        start_time_8=$(date +%s)
        while read -r line; do
                samtools view -L "${line}"_output.bed -h RNAaln_sorted.bam | samtools view -b - > "${line}"_reads.bam &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Bam files for rna reads done" >> ./../../progress.log
        print_elapsed_time $start_time_8 >> ./../../progress.log

        rm *.bed


        start_time_9=$(date +%s)
        while read -r line; do
                bedtools bamtofastq -i "${line}"_reads.bam -fq ./../../rna_fastqs/"${line}".fastq &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Fastq files for rna reads done" >> ./../../progress.log
        print_elapsed_time $start_time_9 >> ./../../progress.log

        cd ..
        cd ..
fi

wait

if [ -n "${DNA_FORWARD}" ]; then

        start_time_8=$(date +%s)

        mkdir dna_Illumina_fastqs
        cd samsANDbams
        mkdir dna_Illumina_alignments
        cd dna_Illumina_alignments

        bwa index ${GENOME}
        bwa mem -t 4 ${GENOME} ${DNA_FORWARD} ${DNA_REVERSE} > DNA_Illumina_aln.sam
        samtools view -Sb ./DNA_Illumina_aln.sam -o ./DNA_Illumina_aln.bam
        samtools sort -o ./DNA_Illumina_aln_sorted.bam ./DNA_Illumina_aln.bam
        samtools index ./DNA_Illumina_aln_sorted.bam ./DNA_Illumina_aln_sorted.bai
        bedtools bamtobed -i ./DNA_Illumina_aln_sorted.bam > output.bed

        max_jobs=12
        job_count=0

        while read -r line; do
                grep "${line}" output.bed > "${line}"_output.bed &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Bed files for DNA Illumina reads done" >> ./../../progress.log
        print_elapsed_time $start_time_8 >> ./../../progress.log

        start_time_9=$(date +%s)
        while read -r line; do
                samtools view -L "${line}"_output.bed -h DNA_Illumina_aln_sorted.bam | samtools view -b - > "${line}"_reads.bam &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Bam files for DNA Illumina reads done" >> ./../../progress.log
        print_elapsed_time $start_time_9 >> ./../../progress.log

        rm *.bed


        start_time_10=$(date +%s)
        while read -r line; do
                bedtools bamtofastq -i "${line}"_reads.bam -fq ./../../DNA_Illumina_fastqs/"${line}".fastq &

                job_count=$((job_count + 1))

                if [ "$job_count" -ge "$max_jobs" ]; then
                        wait
                        job_count=0
                fi

        done < ./../../contig_names.txt

        wait

        echo "Fastq files for DNA Illumina reads done" >> ./../../progress.log
        print_elapsed_time $start_time_10 >> ./../../progress.log

        cd ..
        cd ..
fi

wait




samsANDbams_end=$(date +%s)
samsANDbams_runtime=$((samsANDbams_end - samsANDbams_start))
samsANDbams_runtime_minutes=$((samsANDbams_runtime / 60))
echo -e "bed, bam, and fastq files generated in $samsANDbams_runtime_minutes minutes /n" >> progress.log

sbatch sidrblast.sh
sbatch read_GC.sh
sbatch kmers.sh
sbatch depth.sh
sbatch coverage.sh
sbatch plusANDminus.sh
sbatch contig_GC.sh

#potential future parallelization because all these scripts could run at the same time, but I think it is currently messing it up
