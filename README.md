# Intro

This brief tutorial provides the bare essentials for mapping mutations in pooled NGS data. 
Specifically this pipeline was used to identify a recessive mutation in Tetrahymena.  Files 
that we used are included in this repository.

# Quality control

Illumina paired end reads were checked for quality using fastqc which is avalible below.  The 5' 
and 3' end had lower quality than the rest of the reads.  There was no evidence of primer or adaptor 
contamination. 


   http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Mapping paired end reads to reference genome 

We mapped the reads to against two Tetrahymena thermophila SB210 assemblies:

The Macro nucleous :
   ```
   T_thermophila_June2014_assembly.fasta.txt.zip
   ```
The Micro nucleous : 
   ```
   tetrahymena_thermophila_sb210__mic__2_supercontigs.fasta.zip
   ``` 
**Warning: The SB10 assemblies and annotations continue to improve.  Before you begin your mapping project be sure to check you are working with the most recent assembly and annotation.**

Combine the two provided sequences after unzipping them using cat:
   ```
   cat  T_thermophila_June2014_assembly.fasta.txt > MAC_MIC.fasta
   cat  tetrahymena_thermophila_sb210__mic__2_supercontigs.fasta >> MAC_MIC.fasta
   ```

## Using bowtie2 to map the reads 

Bowtie2 or BWA would both work for this step; we used bowtie2.  Follow the instructions on how to install and use bowtie2.
http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml


Index the reference sequence:

```
bowtie2-build MAC_MIC.fasta MAC_MIC.fasta
```

Align the reads:

```
bowtie2 --met-file B1868 -p 10 --local -x ../ref/MIC.MAC.fasta --n-ceil 1  --trim5 10 --trim3 5 --no-mixed --no-discordant --no-contain --no-overlap  --rg SM:B1868 --rg-id B1868 -1 B1868_ATCACG_L001_R1_001.fastq.gz -2 B1868_ATCACG_L001_R2_001.fastq.gz -S B1868.sam 2> B1868.error
~/tools/bowtie2-2.0.0-beta6/bowtie2 --met-file DisA1 -p 10 --local -x ../ref/MIC.MAC.fasta --n-ceil 1  --trim5 10 --trim3 5 --no-mixed --no-discordant --no-contain --no-overlap  --rg SM:DisA1 --rg-id DisA1 -1 DisA1_TGACCA_L001_R1_001.fastq.gz -2 DisA1_TGACCA_L001_R2_001.fastq.gz -S DisA1.sam 2> DisA1.error
```

