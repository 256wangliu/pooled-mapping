# Intro

This brief tutorial provides the bare essentials for mapping mutations in pooled NGS data. 
Specifically this pipeline was used to identify the recessive *disA* mutation in Tetrahymena.  Files 
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
bowtie2 --met-file B1868 -p 10 --local -x ../ref/MAC_MIC.fasta --n-ceil 1  --trim5 10 --trim3 5 --no-mixed --no-discordant --no-contain --no-overlap  --rg SM:B1868 --rg-id B1868 -1 B1868_ATCACG_L001_R1_001.fastq.gz -2 B1868_ATCACG_L001_R2_001.fastq.gz -S B1868.sam 2> B1868.error
bowtie2 --met-file DisA1 -p 10 --local -x ../ref/MAC_MIC.fasta --n-ceil 1  --trim5 10 --trim3 5 --no-mixed --no-discordant --no-contain --no-overlap  --rg SM:DisA1 --rg-id DisA1 -1 DisA1_TGACCA_L001_R1_001.fastq.gz -2 DisA1_TGACCA_L001_R2_001.fastq.gz -S DisA1.sam 2> DisA1.error
```

   A description of the options used:
   
*  --n-ceil 1  - skip reads with more than one 'N' - these are ambigious bases.  Reads with many ambigious bases should be skipped.
*  --trim5  10     - trim off 10 bases for each read.  This will solely depend on your data.
*  --trim3  5      - trim off 5 bases from each read.  This will solely depend on your data.
*  --no-mixed      - suppress unpaired alignments for paired reads.
*  --no-discordant - suppress discordant alignments for paired reads.
*  --no-contain    - not concordant when one mate alignment contains other.
*  --no-overlap    - not concordant when mates overlap at all.
*  --rg            - sets the sample name. you --rg tag will look like "SM:yourSampleName"
*  --rg-id         - set the read group that corrisoponds to the "SM:yourSampleName".  If you only have one read group, "rg-id" can be the same as "--rg"

## Converting SAM to BAM using Samtools:
```
 ~/tools/samtools-0.1.18/samtools view -bS -f 0x3 B1868.sam > B1868.bam ; rm B1868.sam
 ~/tools/samtools-0.1.18/samtools view -bS -f 0x3 DisA1.sam > DisA1.bam ; rm DisA1.sam
```

The "-bS -f 0x3" flags tell samtools we want to output a bam file "-b" the "S" tells samtools the imput is SAM format.
Only reads that are propoerly paired and mapped properly were kept using "-f 0x3".  For more information on using the "-f" flag look at the samtools documentations and this website:

http://picard.sourceforge.net/explain-flags.html

## Removing optical duplicates.

In this step we remove "optical/pcr duplicates".  These duplicate reads are a product of library preperation.  We use Samtools for this step, but Picard would also work.



removing duplicates:
