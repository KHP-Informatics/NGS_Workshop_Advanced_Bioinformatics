---
title: "Performing variant calling with Freebayes and variant filtering with vcflib"
author: "Alfredo Iacoangeli"
date: "16/03/2020"
---

## IMPORTANT: the paths to files and directories in this workshop are only examples. You will need use your own so please try not to copy and paste the commands but write them yourself matching the correct locations on your computer


Accurate and consistent variant calling requires statistical modelling and is essential for the clinical implementation of NGS. However, many programs are available for calling variants and their concordance varies. Furthermore, variants have different levels of confidence due to differences in data quality. For variants with intermediate confidence levels, it is difficult to separate true variation from artefacts that arise from many factors such as sequencing error, misalignment and inaccurate base quality scores. As a result, the evidence for variant calls requires scrutiny and caution should be used when interpreting positive and negative findings especially for indels which are more error prone. At the end of this exercise you will be able to:

-Use a a state of the art software (Freebayes) to call small variants (SNVs and indels)
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

In its simplest operation, freebayes requires only two inputs: a FASTA reference sequence, and a BAM-format alignment file sorted by reference position. However, the referenced must be indexed. We are going to convert to text format the reference (as required bysamtools faidx), index it with samtools faidx, call variants with freebayes, compress the resultining variant file (VCF) and index the VCF with tabix:

'''
$ zcat ~/ngs_course/dnaseq/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq/data/reference/hg19.fa 

$ samtools faidx ~/ngs_course/dnaseq/data/reference/hg19.fa

$ freebayes --bam ~/ngs_course/dnaseq/data/aligned_data/WES01_chr22m_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq/data/results/WES01_chr22m.vcf

$ bgzip ~/ngs_course/dnaseq/data/results/WES01_chr22m.vcf

$ tabix -p vcf ~/ngs_course/dnaseq/data/results/WES01_chr22m.vcf.gz

'''
