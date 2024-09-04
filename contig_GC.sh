#!/bin/bash
#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_contig_GC.log
#SBATCH --mail-user=vegge003@fiu.edu


contig_GC_start=$(date +%s)

##calculate GC content of contigs, or contig_GC

#titles for columns in the output file
echo -e "contig\tlength\tgc_percent" > ./stats/contig_GC.txt

#loop through each contig in contig_names.txt
while read -r line;
   do
        #grep the line with the contig name, plus the following line | print only the following line |
        #substitute any G or C or g or c to an empty character, counting the number of substitutions and save as the variable gc_number.
        gc_number=$(grep -A 1 "${line}" "${GENOME}" | awk 'NR%2==0 {print}' | awk '{count += gsub(/[gcGC]/, "");} END {print count}')

        #grep the line with the contig name, plus the following line | print only the following line |
        #remove the new line character | count all characters and save as the variable total_bases.
        total_bases=$(grep -A 1 "${line}" "${GENOME}" | awk 'NR%2==0 {print}' | tr -d '\n' | wc -c)

        #pass the variables to awk, divide gc_number by total_bases and multiply by 100
        gc_percent=$(awk -v total="$total_bases" -v value="$gc_number" 'BEGIN { printf((value/total)*100); }')

        #print and append the values for the current contig to the output file
        echo -e "$line\t$total_bases\t$gc_percent" >> ./stats/contig_GC.txt
   done < contig_names.txt

cd ..

contig_GC_end=$(date +%s)
contig_GC_runtime=$((contig_GC_end - contig_GC_start))
contig_GC_runtime_minutes=$((contig_GC_runtime / 60))
echo "contig_GC stats completed in $contig_GC_runtime_minutes minutes" >> progress.log
