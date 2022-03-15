---
title: "NGS Variant Annotation"
author: "Alfredo Iacoangeli"
date: "15/03/2021"
---

## Variant Annotation and Prioritisation

At this point in the analysis of whole exome data from patient WES01, we have assessed the quality of raw sequence data, removed low quality reads, aligned sequence data to the reference genome, removed poorly mapped reads and duplicate reads, assessed the alignment process, called variants and evaluated variant calling. The aim of this practical is to complete the analysis of WES01 by annotating the variants and using a filtering strategy to generate a short-list of potentially pathogenic variants. At the end of this exercise you will be able to:

Use ANNOVAR to annotate variants with respect to genes, databases of normal variation, and predictors of pathogenicity
Apply basic variant filters to generate a list of candidate genes and likely causative variants
Use phenotype information to improve annotation, filtering and prioritisation of likely causative variants

### Annotating variants
The result of variant calling is a variant call file (VCF) which describes the location and data quality of each variant. However, the initial VCF file does not provide any information about the functional relevance of variants and their potential contribution to disease. To gain these insights we will use ANNOVAR (Wang et al 2010) to annotate each variant in the VCF file with respect to their location in relation to genes and coding sequences (exon, intron or intergenic), whether they change the coding sequence and if so how (missense, stop gain, synonymous, frameshift, amino acid consequence etc). In addition, we will cross reference the variants with databases of known variation (1000 genomes, dbSNP, Exome Sequencing Project and COSMIC) and predictions of functional significance (avsift and conservation scores). At this stage, it is important to be aware that the annotation result will vary according to the choice of annotation software (alternative software: SnpEff Cingolani et al 2012 and Variant Effect Predictor VEP McLaren et al 2010) and definition of transcripts (Ensemble Flicek et al 2012, RefSeq Pruitt et al 2012, Gencode Harrow et al 2012). This variability occurs because a single variant may have different effects in different transcripts (isoforms of the gene) and in some cases may affect different genes (genes can overlap, one on forward strand the other on the reverse strand) each with multiple transcripts. It is also possible to interpret a single variant as having multiple effects on a single transcript (Figure 1). The annotation software must therefore apply a set of rules to determine which variant takes precedence see here for ANNOVAR precedence rules but these rules vary between programs so that they generate different results. For example, the overall agreement for Loss of Function (LoF) variants annotated by ANNOVAR and VEP is just 64% (McCarthy et al 2014 : Choice of transcripts and software has a large effect on variant annotation).
The choice of which transcript set to use for annotation has an even bigger effect on the annotation results. For example, the overlap between LoF variants is only 44% when the same software is used to annotate against transcripts from Ensembl or RefSeq (McCarthy et al 2014). This is because RefSeq transcripts are based on a collection of non-redundant mRNA models that have strong support and are manually curated. As a result, the RefSeq database is highly accurate but does not contain all possible transcripts or gene models (n=41,501 transcripts from RefSeq release 57 that were used by ANNOVAR). In comparison, Ensembl provides a less curated but more comprehensive list of transcripts (n=115,901 transcripts from Ensembl version 69 that were used by ANNOVAR) based on information from the Consensus Coding Sequence (CCDS, Pruitt et al 2009), Human and Vertebrate Analysis and Annotation (HAVANA), Vertebrate Genome Annotation (Vega, Ashurst et al 2005), ENCODE (The ENCODE Project Consortium 2012) and GENCODE (Searle et al 2010).

ANNOVAR is a rapid, efficient tool to annotate functional consequences of genetic variation from high-throughput sequencing data. wANNOVAR provides easy and intuitive web-based access to the most popular functionalities of the ANNOVAR software

This tool will annotate variants using specified gene annotations, regions, and filtering databases. Input is a VCF dataset, and output is a table of annotations for each variant in the VCF dataset or a VCF dataset with the annotations in INFO fields.

ANNOVAR Website: http://annovar.openbioinformatics.org
Paper: http://nar.oxfordjournals.org/content/38/16/e164

### Download Annovar

Annovar is not opensource and we could not install it with the conda command. However, it is free for academic use.
Please register on the Annovar website and request a copy for this workshop https://www.openbioinformatics.org/annovar/annovar_download_form.php
After registering you will recive an email with a link to download the source code of Annovar on your personal computer. Please click on the link in the email, download the Annovar source code and upload on your Openstack instance using FileZilla (https://filezilla-project.org).
We have used FileZilla before in a previous part of the practical (https://github.com/KHP-Informatics/NGS_Workshop_Advanced_Bioinformatics/tree/master/assessing_quality), please see the "assessing_quality" workshop if you do not remember how to use FileZilla to download/upload files to/from your OpenStack instance.

After uploading the annovar.latest.tar.gz file onto OpenStack, please unpack it and set annovar up:

```
$ tar -zxvf annovar.latest.tar.gz
```


### Download Annovar db

Before using Annovar you need to download the databases it uses for the annotation. The following commands will download some nasic ones: 

```
$ cd annovar
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
$ ./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
```

#### If you get a "permission denied" error message when you try to download the DBs, this might be due to annotate_variation.pl not to be set as an executable file. Please fix this but running "chmod +x annotate_variation.pl" in your temrinal

### VCF to Annovar input format

```bash
$ ./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.vcf.gz > ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.avinput
```
```bash
NOTICE: Finished reading N lines from VCF file
NOTICE: A total of N locus in VCF file passed QC threshold, representing N SNPs (N transitions and N transversions) and N indels/substitutions
NOTICE: Finished writing N SNP genotypes (N transitions and N transversions) and N indels/substitutions for 1 sample
```

### Run Annovar table function

* csv output
```bash
$ ./table_annovar.pl ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22.avinput humandb/ -buildver hg19  \ 
   -out ~/ngs_course/dnaseq/results/WES01_chr22m_filtered_chr22 -remove   \ 
      -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout
```

The output will be in csv format, please download it via FileZilla and open it with Office Excel or any other speadsheet software. 


## Exercise

Take a moment to update the README for the dnaseq folder (hint: use vim, or nano or any text editor of your choice to create the file). Give a short update of the project and brief descriptions of the types of file you have generated within each of the sub-directories. Please take note of the current size of the project. The total storage available on your virtual machine if 40Gigabytes

* Reference by:
    * [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
    * [Using VAAST to Identify an X-Linked Disorder Resulting in Lethality in Male Infants Due to N-Terminal Acetyltransferase Deficiency](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135802/)
    * [ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data](http://www.ncbi.nlm.nih.gov/pubmed/20601685)
    * [Genomic variant annotation and prioritization with ANNOVAR and wANNOVAR](http://www.nature.com/nprot/journal/v10/n10/full/nprot.2015.105.html)
    * [NAA10 / rs387906701](http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=rs387906701)
    * [SNPedia](http://snpedia.com/index.php/Rs387906701)
    * [NAA10 mutation causing a novel intellectual disability syndrome with Long QT due to N-terminal acetyltransferase impairment](http://www.nature.com/articles/srep16022)
    * [OMIM - Ogden Syndrome](http://www.omim.org/entry/300855)

