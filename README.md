# Intro

This brief tutorial provides the bare essentials for mapping mutations in pooled NGS data. 
Specifically this pipeline was used to identify the recessive *DisA* mutation in *Tetrahymena*.  Please have 
a look at the paper for details about the experimental design and nature of  *DisA*.
There are many ways to map, call, and identify mutations.  This documents ONLY describes one successful approach. 
Diverge from this document as it becomes outdated or at your fancy.  


# Quality control

Fastqc was used to quality control the Illumina paired-end reads.  The 5' 
and 3' end had lower quality than the rest of the reads.  There was no evidence of primer or adaptor 
contamination. 


   http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Mapping paired end reads to reference genome 

We mapped the reads to against two *Tetrahymena thermophila* SB210 assemblies:

The Macro nucleous :
   ```
   T_thermophila_June2014_assembly.fasta.txt.zip
   ```
The Micro nucleous : 
   ```
   tetrahymena_thermophila_sb210__mic__2_supercontigs.fasta.zip
   ``` 
**Warning: The SB10 assemblies and annotations continue to improve.  Before you begin your mapping project be sure to check that you are working with the most recent assembly and annotation.**

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

**Warning: It is important to understand the rg and rg-id settings.  Failures to properly specify either flag will result in a massive headache downstream!**

## Converting SAM to BAM using Samtools:
```
 ~/tools/samtools-0.1.18/samtools view -bS -f 0x3 B1868.sam > B1868.bam ; rm B1868.sam
 ~/tools/samtools-0.1.18/samtools view -bS -f 0x3 DisA1.sam > DisA1.bam ; rm DisA1.sam
```

The "-bS -f 0x3" flag tell samtools to output a bam file "-b" and the input is SAM format "S".
The "-f 0x3" flag only keeps reads that are paired and properly mapped.  For more information on using the "-f" flag look at the samtools documentation and this website:

http://picard.sourceforge.net/explain-flags.html

## Sorting the BAM files by start coordinates
```
~/tools/samtools-0.1.18/samtools sort  B1868.bam  B1868.sort
~/tools/samtools-0.1.18/samtools sort  DisA1.bam  DisA1.sort
```

## Removing optical duplicates

In this step we remove "optical/pcr duplicates".  These duplicate reads are a product of library preperation.  We use Samtools for this step, but Picard would also work.

removing duplicates:

```
 ~/tools/samtools-0.1.18/samtools rmdup B1868.sort.bam B1868.sort.rmdup.bam
 ~/tools/samtools-0.1.18/samtools rmdup DisA1.sort.bam DisA1.sort.rmdup.bam
 ```
 
## Indexing the BAMs

 ```
 ~/tools/samtools-0.1.18/samtools index B1868.sort.rmdup.bam
 ~/tools/samtools-0.1.18/samtools index DisA1.sort.rmdup.bam
 ```

## Indel realignment 

Indel re-alignment is a two step process that improves the alignments in regions with many INDELs.  In the first step GATK identifies intervals that need to be re-aligned and in the second step they are realigned.  Much more information can be found at the GATK website:
https://www.broadinstitute.org/gatk/

creating target intervals:
```
java -Xmx10g -jar /usr/local/GenomeAnalysisTK-2.4-7-g5e89f01/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../ref/MIC.MAC.fasta -I  B1868.sort.rmdup.bam -o B1868.sort.rmdup.intervals
java -Xmx10g -jar /usr/local/GenomeAnalysisTK-2.4-7-g5e89f01/GenomeAnalysisTK.jar -T RealignerTargetCreator -R ../ref/MIC.MAC.fasta -I  DisA1.sort.rmdup.bam -o DisA1.sort.rmdup.intervals
```

running the realignment:
```
java -Xmx20g -jar /usr/local/GenomeAnalysisTK-2.4-7-g5e89f01/GenomeAnalysisTK.jar -T IndelRealigner  -R ../ref/MIC.MAC.fasta -I  B1868.sort.rmdup.bam -targetIntervals B1868.sort.rmdup.intervals -o B1868.sort.rmdup.realn.bam
java -Xmx20g -jar /usr/local/GenomeAnalysisTK-2.4-7-g5e89f01/GenomeAnalysisTK.jar -T IndelRealigner  -R ../ref/MIC.MAC.fasta -I  DisA1.sort.rmdup.bam -targetIntervals DisA1.sort.rmdup.intervals -o DisA1.sort.rmdup.realn.bam
```

## Calling Variants

For variant calling we used SNVer pooled.  Again there are several other programs that would work for variant calling.  It very much depends on what kind of sample you are working on.  For example if you are trying to identify low frequency variants the tool *LoFreq* might be right for you.  SNVer pooled produces a file for Single Nucleotide Variants SNVs and INDELs.  

The config file "snver.file.txt":
```
#name	no.haploids	no.samples	mq	bq
B1868.sort.rmdup.realn.bam	2	1	20	20
DisA1.sort.rmdup.realn.bam	32	16	20	20

```

This information is use by SNVER statistical model.  Notice the *DisA1* sample acturally consists of 16 clonal samples (See paper for more information).

The command line to call variants:

```
java -Xmx50g -jar ~/tools/SNVerPool.jar -i ../raw_reads/ -t 0 -a 0 -u 2 -r ../ref/MIC.MAC.fasta -c snver.file.txt -o teth
```


## Finding canidates by filtering samples based on allele frequency

We knew *DisA* was recessive so we setup simple heuristic filters.  This was done with a provided script "filter-recessive.pl."

getting help:
```
perl filter-recessive.pl --help
```
 usage statement:
 
 ```
 Synopsis:

cat my-snver-output.vcf | filter-recessive.pl

Description:

This script filters a SNVER pooled file.

WARNING: the recessive mutant must be the first genotype column followed by the
reference strain.  Any other ordering will not work.  Additional columns will be ignored

Output:

    seqid                      : The contig / scaffold / chromosome name
    position                   : The position of the mutation
    reference allele frequency : The frequency of the non-reference mutation in the control pool
    mutant alllele frequency   : The frequency of the non-reference mutation in the mutant pool
    delta af                   : The allele frequency difference between the mutant and control pool
    fail - mutant in control   : We do not expect to see the mutation in the reference pool [0 = pass, 1 = fail]
    fail - mutant af too low   : Under the recessive model we expect the allele frequency to be 1.0, however
                                 sequencing errors could result in lower allele frequencies so the cutoff is 0.75 [0 = pass, 1 = fail]
    fail - depth               : If the control or mutant pool has a depth below 5 [0 = pass, 1 = fail]
    fail flags                 : The last three columns concatenated

Contact:

   Problems or questions: zev.kronenberg@gmail.com
 ```
