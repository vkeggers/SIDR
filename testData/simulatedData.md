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
module load samtools-1.15.1-gcc-8.2.0
cat *.sam > HIFIelegans.sam   #concatenate all the sam files into one
samtools view -b ./HIFIelegans.sam -o ./HIFIelegans.bam
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

ccs HIFIelegans.bam HIFIelegans.fastq.gz
```

**Step 4: repeat** pbsim3, samtools, ccs (steps 1-3) for ecoli, changing the genome from caenorhabditis_elegans.PRJNA13758.WBPS18.genomic.fa
to GCF_000008865.2_ASM886v2_genomic.fna and output file names from HIFIelegans to HIFIecoli

**Step 5: unzip and concatenate the fastq files together:**

```
gunzip HIFIelegans.fastq.gz
gunzip HIFIecoli.fastq.gz
cat HIFIelegans.fastq HIFIecoli.fastq > HIFItestData.fastq
```

JUST KIDDING :<

When you go to assemble this the assembly is only 20M with 2 contigs??? Let's try using less ecoli sequences by using seqtk to take 100 random sequences:

```
module load seqtk-1.3
seqtk sample -s100 HIFIecoli.fastq 100 > subsetHIFIecoli.fastq
```

</details>


<details>
<summary>ONT</summary>
  
</details>


<details>
<summary>Paired-end Illumina</summary>
  
</details>


### RNA data acquisition and generating fastq:

<details>
<summary>Data</summary>

```
#from wormbase
wget https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS18/species/caenorhabditis_elegans/PRJNA13758/caenorhabditis_elegans.PRJNA13758.WBPS18.CDS_transcripts.fa.gz
gunzip caenorhabditis_elegans.PRJNA13758.WBPS18.CDS_transcripts.fa.gz
```

As for the _E.coli_ the same GCF_000008865.2_ASM886v2_genomic.fna assembly file was used. 

  
</details>




<details>
<summary>fastq</summary>
  
</details>



### Statistics for raw data calculated with fastqc:
