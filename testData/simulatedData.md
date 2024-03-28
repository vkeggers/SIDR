### Reference assemblies for _Caenorhabditis_ and _E.coli_ used in generating test data:

1. caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa.gz

2. GCF_000008865.2_ASM886v2_genomic.fna

### Data acquisition:

<details>
<summary><i>C.elegans</i></summary>

```
#from wormbase
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa.gz
gunzip caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa.gz
```
  
</details>



<details>
<summary><i>E.coli</i></summary>

```
#from NCBI
module load sratoolkit-3.0.0
fastq-dump GCF_000008865.2
unzip ncbi_dataset\ \(1\).zip
cd ncbi_dataset/data/GCF_000008865.2
cp GCF_000008865.2_ASM886v2_genomic.fna ./../../../.
rm -r ncbi_dataset
```
  
</details>




### How fastq files were generated from reference assemblies:

<details>
<summary>HIFI</summary>
  
The HIFI fastq files were made with the software [pbsim3](https://github.com/yukiteruono/pbsim3) in multipass mode.
This creates a sam file which must be converted into a bam using the software [samtools](https://www.htslib.org/).
This bam is then input into pacbio's [ccs](https://ccs.how/) software, which was installed with the bioconda package pbccs.

**Step 1: pbsim3**

```
#!/bin/bash
#SBATCH --account account_name
#SBATCH --qos partition_name
#SBATCH --partition partition_name
#SBATCH --output=out_%pbsim.log

module load pbsim3-3.0.4

pbsim --strategy wgs \
      --method qshmm \
      --qshmm QSHMM-RSII.model \
      --depth 60 \
      --genome /your/path/to/caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa \
      --pass-num 10
```

This will generate a sam file for each chromosome (for elegans it is 7, 6 plus mitochondria)


**Step 2: samtools**

```
#!/bin/bash
#SBATCH --account account_name
#SBATCH --qos partition_name
#SBATCH --partition partition_name
#SBATCH --output=out_%_sam.log

module load samtools-1.9-gcc-8.2.0-o53igvd

ls *.sam > list.txt

while read -r file; do
    samtools view -b "$file" -o "$file".bam
done < list.txt
```

**Step 3: ccs**

```
#!/bin/bash
#SBATCH --account account_name
#SBATCH --qos partition_name
#SBATCH --partition partition_name
#SBATCH --output=out_%ccs.log

module load mamba/23.1.0-4
source activate pbccs

rm list.txt

ls *.bam > list.txt

while read -r file; do
    ccs "$file" "$file".fastq.gz
done < list.txt
```

**Step 4: unzip and concatenate the fastq files together:**

```
gunzip *.fastq
cat *.fastq > HIFIelegans.fastq
```

**Step 5: repeat** pbsim3, samtools, ccs (steps 1-3) for ecoli, changing the genome from caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa
to GCF_000008865.2_ASM886v2_genomic.fna and output file name from HIFIelegans.fastq to HIFIecoli.fastq

**Step 6: concatenate the fastq files together:** 

```
cat HIFIelegans.fastq HIFIecoli.fastq > HIFItestData.fastq
```


</details>


<details>
<summary>ONT</summary>

The HIFI fastq files were made with the software [pbsim3](https://github.com/yukiteruono/pbsim3) with the high quality ONT model, generating reads with at least 90% accuracy. 

```
#!/bin/bash
#SBATCH --account account_name
#SBATCH --qos partition_name
#SBATCH --partition partition_name
#SBATCH --output=out_%pbsim.log

module load pbsim3-3.0.4

pbsim --strategy wgs \
      --method qshmm \
      --qshmm QSHMM-ONT-HQ.model \
      --depth 60 \
      --genome ./caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa
```

repeat w/ e.coli genome and then concatenate all .fastq files together under ONTtestData.fastq

</details>


<details>
<summary>Paired-end Illumina</summary>

All Illumina data was generated using the software [ART](https://manpages.debian.org/testing/art-nextgen-simulation-tools/art_illumina.1.en.html), which was installed with the bioconda package.

```
#!/bin/bash

#SBATCH --account account_name
#SBATCH --qos partition_name
#SBATCH --partition partition_name
#SBATCH --output=out_art.log
#SBATCH --mail-user=username@email.com   #use your own email
#SBATCH --mail-type=ALL
#SBATCH -n 8

module load mamba/23.1.0-4
source activate art

art_illumina -sam -i caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa -l 150 -p -nf 0 -f 60 -m 200 -s 10 -ss HS25 -o Illumina_elegans
```

repeat with ecoli and concatenate the Illumina_elegans1 with Illumina_ecoli1 and Illumina_elegans2 with Illumina_ecoli2
  
</details>

<details>
<summary>RNA Illumina</summary>

<details>
<summary>Data acquisition</summary>
  
```
#from wormbase
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS18.CDS_transcripts.fa.gz
gunzip caenorhabditis_elegans.PRJNA13758.WBPS18.CDS_transcripts.fa.gz
```

As for the _E.coli_ the same GCF_000008865.2_ASM886v2_genomic.fna assembly file was used. 

</details>



<details>
<summary>fastq generation</summary>

Repeat of paired-end Illumina but this time usingcaenorhabditis_elegans.PRJNA13758.WBPS18.CDS_transcripts.fa instead of caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa

Again concatenate this with the ecoli data.
  
</details>



</details>







  

  
</details>



### Statistics for raw data calculated with fastqc:

<details>
<summary>TABLE 1: files and files sizes</summary>
  
| file | size | number of sequences |
|------|------|---------------------|
| caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa | 98M | 7 |
| caenorhabditis_elegans.PRJNA13758.WBPS18.CDS_transcripts.fa | 41M | 28,558 |
| GCF_000008865.2_ASM886v2_genomic.fna | 5.5M | 3 |
| HIFItestData.fastq | NA | NA |
| ONTtestData.fastq | 12G | 689,630 |
| IlluminaTestData_1.fastq | 6.3G | 21,176,070 |
| IlluminaTestData_2.fastq | 6.3G | 21,176,070 |
| IlluminaRNAtestData_1.fastq | NA | NA |
| IlluminaRNAtestData_2.fastq | NA | NA |

</details>

<details>
<summary>TABLE 2: Quality score, read lengths, and GC%</summary>

| file | mean Phred Quality Score | shortest read | longest read | mean read length | GC% |
|------|--------------------------|---------------|--------------|------------------|-----|
| HIFItestData.fastq | NA | NA | NA | NA | NA |
| ONTtestData.fastq | 14 | 96 | 78,823 | 8000 | 36% |
| IlluminaTestData_1.fastq | 36 | 150 | 150 | 150 | 36% |
| IlluminaTestData_2.fastq | 36 | 150 | 150 | 150 | 36% |
| IlluminaRNAtestData_1.fastq | NA | NA | NA | NA | NA |
| IlluminaRNAtestData_2.fastq | NA | NA | NA | NA | NA |

  
</details>


### Can we make contaminated assemblies?

<details>
<summary>hifiasm</summary>

</details>


<details>
<summary>flye with canu correct</summary>

A genome assembly of ONT reads was created with the ONTtestData.fastq file using flye with the canu correct module as done in https://github.com/vkeggers/DeNovo-Nematode-Pipeline/blob/main/DeNovo-Nematode-Pipeline.md

This assembly was blasted against the nt database and 2 of the 15 contigs hit to _e.coli_, showing that the contaminant reads were incorporated into the genome assembly. This assembly was polished with the IlluminaTestData_1.fastq and IlluminaTestData_2.fastq using pilon, again following  https://github.com/vkeggers/DeNovo-Nematode-Pipeline/blob/main/DeNovo-Nematode-Pipeline.md

The polished genome assembly was again blasted against nt.


</details>



<details>
<summary>nextdenovo</summary>
  
</details>
