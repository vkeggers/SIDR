#!/bin/bash

#SBATCH --account iacc_jfierst
#SBATCH --qos highmem1
#SBATCH --partition highmem1
#SBATCH --output=out_blastnt.log
#SBATCH --mail-user=your@email.com
#SBATCH --mail-type=ALL

module load blast-plus-2.7.1-gcc-8.2.0-fm5yf3k


#makeblastdb -in /home/data/jfierst/veggers/nt_db/nt.fa -dbtype nucl -out nt

GENOME=/home/data/jfierst/veggers/Oscheius/DF5033/DF5033_hifiAssembly/DF5033ONTpb.asm.bp.p_ctg.fa
NT=/home/data/jfierst/veggers/nt_db/nt

blastn -query ${GENOME} \
    -db ${NT} \
   # -outfmt 6 \
    -culling_limit 5 \
    -evalue 1e-25 \
    -out SIDRblast.txt


################################
################################
########FIX_BLAST_OUTPUT########
################################
################################

#blast results are by default sorted with the most significant e-value occurring first
#get the first blast results for each contig into a list
awk '/Query=/{print; flag=1; next} flag && /^>/{print; flag=0} flag && /No hits found/{print; flag=0}' out_blastnt.log > test1.txt
awk 'NR%2==1{col1=$0} NR%2==0{print col1, $0}' test1.txt > test2.txt
sed -i 's/PREDICTED: //' test2.txt
awk '{print $2, $5}' test2.txt > test3.txt
sed -i 's/hits found/No_hits_found/g' test3.txt
sed -i 's/ /\t/' test3.txt
echo -e "contig\tTopHit" > TopHitBLAST.txt
sort -k1 test3.txt >> TopHitBLAST.txt


#get the most frequent blast results into a list
awk '/Query=/{print; flag=1; next} flag && /^>/{print; next} flag && /No hits found/{print; flag=0}' testData_flye_blast.txt > test1.txt
sed 's/PREDICTED: //' test1.txt > test2.txt
awk '/Query=/{print; flag=1; next} flag && /^>/{print $1,$3; next} flag && /No hits found/{print; flag=0}' test2.txt > test3.txt
sed 's/\*\*\*\*\* No hits found \*\*\*\*\*/> No_hits_found/g' test3.txt > test4.txt
awk '/Query=/{query=$0; next} $0 ~ /^>/ {result[query][$0]++} END {for (q in result) {max_count=0; most_frequent=""; for (r in result[q]) {if (result[q][r] > max_count) {max_count=result[q][r]; most_frequent=r}} print q, most_frequent}}' test4.txt > test5.txt
awk '{print $2,$4}' test5.txt > test6.txt
sed -i 's/ /\t/' test6.txt
echo -e "contig\tMostHits" > MostFrequentBLAST.txt
sort -k1 test6.txt >> MostFrequentBLAST.txt

#remove all temporary files
rm test*.txt
