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
BLAST_OUTPUT=/home/data/jfierst/veggers/testData_flye_blast.txt

blastn -query ${GENOME} \
    -db ${NT} \
   # -outfmt 6 \
    -culling_limit 5 \
    -evalue 1e-25 \
    -out BLASToutput.txt


################################
################################
########FIX_BLAST_OUTPUT########
################################
################################

#blast results are by default sorted with the most significant e-value occurring first
#get the first blast results for each contig into a list
awk '/Query=/{print; flag=1; next} flag && /^>/{print; flag=0} flag && /No hits found/{print; flag=0}' ${BLAST_OUTPUT} > temp1.txt
awk 'NR%2==1{col1=$0} NR%2==0{print col1, $0}' temp1.txt > temp2.txt
sed -i 's/PREDICTED: //' temp2.txt
sed -i 's/\*\*\*\*\* No hits found \*\*\*\*\*/> No No_hits_found/g' temp2.txt 
awk '{print $2, $5}' temp2.txt > temp3.txt
sed -i 's/ /\t/' temp3.txt
echo -e "contig\tTopHit" > TopHitBLAST.txt
sort -t_ -k2n temp3.txt >> TopHitBLAST.txt


#get the most frequent blast results into a list
awk '/Query=/{print; flag=1; next} flag && /^>/{print; next} flag && /No hits found/{print; flag=0}' ${BLAST_OUTPUT} > temp1.txt
sed 's/PREDICTED: //' temp1.txt > temp2.txt
awk '/Query=/{print; flag=1; next} flag && /^>/{print $1,$3; next} flag && /No hits found/{print; flag=0}' temp2.txt > temp3.txt
sed 's/\*\*\*\*\* No hits found \*\*\*\*\*/> No_hits_found/g' temp3.txt > temp4.txt
awk '/Query=/{query=$0; next} $0 ~ /^>/ {result[query][$0]++} END {for (q in result) {max_count=0; most_frequent=""; for (r in result[q]) {if (result[q][r] > max_count) {max_count=result[q][r]; most_frequent=r}} print q, most_frequent}}' temp4.txt > temp5.txt
awk '{print $2,$4}' temp5.txt > temp6.txt
sed -i 's/ /\t/' temp6.txt
echo -e "contig\tMostHits" > MostFrequentBLAST.txt
sort -t_ -k2n temp6.txt >> MostFrequentBLAST.txt

#remove all temporary files
rm temp*.txt

#join the two blast outputs together
join --check-order --header -t$'\t' TopHitBLAST.txt MostFrequentBLAST.txt > SIDRblast.txt

#I'm not sure that this is worth keeping because there were no instances when the most hit was Oscheius but the top hit wasn't
