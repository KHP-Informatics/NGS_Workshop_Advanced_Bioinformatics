---
title: "Performing variant calling with Freebayes and variant filtering with vcflib"
author: "Alfredo Iacoangeli"
date: "15/03/2021"
---

## IMPORTANT: the paths to files and directories in this workshop are only examples. You will need use your own so please try not to copy and paste the commands but write them yourself matching the correct locations on your computer


Accurate and consistent variant calling requires statistical modelling and is essential for the clinical implementation of NGS. However, many programs are available for calling variants and their concordance varies. Furthermore, variants have different levels of confidence due to differences in data quality. For variants with intermediate confidence levels, it is difficult to separate true variation from artefacts that arise from many factors such as sequencing error, misalignment and inaccurate base quality scores. As a result, the evidence for variant calls requires scrutiny and caution should be used when interpreting positive and negative findings especially for indels which are more error prone. At the end of this exercise you will be able to:

-Use a state of the art software (Freebayes) to call small variants (SNVs and indels)
-Describe the contents of variant call files
-Generate and interpret basic variant quality control parameters
-Use quality control filters to exclude or flag variants with low confidence
-See: https://wiki.galaxyproject.org/Learn/GalaxyNGS101#Finding_variants

## Variant calling with Freebayes


FreeBayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment. See https://github.com/ekg/freebayes for details on FreeBayes.

This should be your directory tree:

```
dnaseq
	├── data
	│   ├── reference
	│   │   └── hg19.fa.gz
	│   │   └── hg19.fa.gz.amb
	│   │   └── etc
 	|   ├── untrimmed_fastq
	│   │   └── WES01_chr22m_R1.fastq.gz
	│   │   └── etc
	│   ├── trimmed_fastq
	│   │   ├── WES01_chr22m_trimmed_R_1P.fastq
	│   │   ├── etc     
	│   └── aligned_data
	│       └── WES01_chr22m_sorted_filtered.bam
	│       └── etc
	├── meta
	├── results
	└── logs
```

FreeBayes uses short-read alignments (BAM files with Phred+33 encoded quality scores, now standard) for any number of individuals from a population and a reference genome (in FASTA format) to determine the most-likely combination of genotypes for the population at each position in the reference. It reports positions which it finds putatively polymorphic in variant call file (VCF) format. It can also use an input set of variants (VCF) as a source of prior information, and a copy number variant map (BED) to define non-uniform ploidy variation across the samples under analysis.

In its simplest operation, freebayes requires only two inputs: a FASTA reference sequence, and a BAM-format alignment file sorted by reference position. However, the referenced must be indexed. We are going to convert to text format the reference (as required bysamtools faidx), index it with samtools faidx, call variants with Freebayes, compress the resulting variant file (VCF) and index the VCF with tabix:

```
$ zcat ~/ngs_course/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa 

$ samtools faidx ~/ngs_course/dnaseq/data/reference/hg19.fa

$ freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/WES01_chr22m_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq/results/WES01_chr22m.vcf

$ bgzip ~/ngs_course/dnaseq/results/WES01_chr22m.vcf

$ tabix -p vcf ~/ngs_course/dnaseq/results/WES01_chr22m.vcf.gz

```


## Filtering the VCF

As we run Freebayes with default parameters, the resulting VCF contains a large number of "bad" calls. The Freebayes Information Fields we will use for filtering are:

QUAL=probability that there is a polymorphism at the loci described by the record: 1 - P(locus is homozygous given the data).
AO=Alternate allele observations, with partial observations recorded fractionally
SAF=Number of alternate observations on the forward strand
SAR=Number of alternate observations on the reverse strand
RPL=Reads Placed Left: number of reads supporting the alternate balanced to the left (5’) of -the alternate allele
RPR=Reads Placed Right: number of reads supporting the alternate balanced to the right (3’) of the alternate allele
What calls do we “know” are poor (w.r.t. freebayes VCF annotations)?
- low-quality calls (QUAL < N)
+ also, maybe QUAL is high but QUAL / AO is low
- loci with low read depth (DP < N)
- alleles that are only seen on one strand
+ remove by “SAF > 0 & SAR > 0”
- alleles that are only observed by reads placed to the left or right
+ remove by “RPL > 0 & RPR > 0”

We will apply the following suggested freebayes hard filter for human diploid sequencing and use vcffilter to apply it:

QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1

QUAL > 1: removes horrible sites
QUAL / AO > 10 : additional contribution of each obs should be 10 log units (~ Q10 per read)
SAF > 0 & SAR > 0 : reads on both strands
RPR > 1 & RPL > 1 : at least two reads “balanced” to each side of the site

```
$ conda install vcflib #if you did not install vcflib before. vcffilter is part of the vcflib suite

$ vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" \ 			
	~/ngs_course/dnaseq/results/WES01_chr22m.vcf.gz > ~/ngs_course/dnaseq/results/WES01_chr22m_filtered.vcf
```

The bed file chr22.genes.b37.bed describes the exome sequences and genes that have been targeted in your trial data. Using bedtools we can filter the vcf file for the regions in chr22.genes.b37.bed: 

```

$ bedtools intersect -header -wa -a ~/ngs_course/dnaseq/results/WES01_chr22m_filtered.vcf -b ../chr22.genes.hg19.bed \
 	> ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.vcf

$ bgzip ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.vcf

$ tabix -p vcf ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.vcf.gz

```


### Exercise

Take a moment to update the README for the dnaseq folder (hint: use vim, or nano or any text editor of your choice to create/edit the file). Give a short update of the project and brief descriptions of the types of file you have generated within each of the sub-directories. Please take note of the current size of the project. The total storage available on your virtual machine if 40Gigabytes
