# Intro

This brief tutorial provides the bare essentials for mapping mutations in pooled NGS data. 
Specifically this pipeline was used to identify a recessive mutation in Tetrahymena.  Files 
that we used are included in this repository.

# Quality control

Illumina paired end reads were check for quality using fastqc which is avalible below.  The 5' 
and 3' end had lower quality than the rest of the reads.  There was no evidence of primer or adaptor 
contamination. 


   http://www.bioinformatics.babraham.ac.uk/projects/fastqc/

# Mapping paired end reads to 

We mapped the reads to against two Tetrahymena thermophila SB210 assemblies:

The Macro nucleous :

The Micro nucleous : 

**Warning: The SB10 assemblies and annotations continue to improve.  Before you begin your mapping project be sure to check you are working with the most recent assembly and annotation.**
