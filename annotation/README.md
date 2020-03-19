# Advanced Bioinformatics module workshop

### NGS Variant Annotation
---
* Date: 16/03/2020
* Place: KCL




## Variant Annotation and Prioritisation

At this point in the analysis of whole exome data from patient WES01, we have assessed the quality of raw sequence data, removed low quality reads, aligned sequence data to the reference genome, removed poorly mapped reads and duplicate reads, assessed the alignment process, called variants and evaluated variant calling. The aim of this practical is to complete the analysis of WES01 by annotating the variants and using a filtering strategy to generate a short-list of potentially pathogenic variants. At the end of this exercise you will be able to:

Use ANNOVAR to annotate variants with respect to genes, databases of normal variation, and predictors of pathogenicity
Apply basic variant filters to generate a list of candidate genes and likely causative variants
Use phenotype information to improve annotation, filtering and prioritisation of likely causative variants

### Annotating variants
The result of variant calling is a variant call file (VCF) which describes the location and data quality of each variant. However, the initial VCF file does not provide any information about the functional relevance of variants and their potential contribution to disease. To gain these insights we will use ANNOVAR (Wang et al 2010) to annotate each variant in the VCF file with respect to their location in relation to genes and coding sequences (exon, intron or intergenic), whether they change the coding sequence and if so how (missense, stop gain, synonymous, frameshift, amino acid consequence etc). In addition, we will cross reference the variants with databases of known variation (1000 genomes, dbSNP, Exome Sequencing Project and COSMIC) and predictions of functional significance (avsift and conservation scores). At this stage, it is important to be aware that the annotation result will vary according to the choice of annotation software (alternative software: SnpEff Cingolani et al 2012 and Variant Effect Predictor VEP McLaren et al 2010) and definition of transcripts (Ensemble Flicek et al 2012, RefSeq Pruitt et al 2012, Gencode Harrow et al 2012). This variability occurs because a single variant may have different effects in different transcripts (isoforms of the gene) and in some cases may effect different genes (genes can overlap, one on forward strand the other on the reverse strand) each with multiple transcripts. It is also possible to interpret a single variant as having multiple effects on a single transcript (Figure 1). The annotation software must therefore apply a set of rules to determine which variant takes precedence see here for ANNOVAR precedence rules but these rules vary between programs so that they generate different results. For example, the overall agreement for Loss of Function (LoF) variants annotated by ANNOVAR and VEP is just 64% (McCarthy et al 2014 : Choice of transcripts and software has a large effect on variant annotation).
The choice of which transcript set to use for annotation has an even bigger effect on the annotation results. For example, the overlap between LoF variants is only 44% when the same software is used to annotate against transcripts from Ensembl or RefSeq (McCarthy et al 2014). This is because RefSeq transcripts are based on a collection of non-redundant mRNA models that have strong support and are manually curated. As a result, the RefSeq database is highly accurate but does not contain all possible transcripts or gene models (n=41,501 transcripts from RefSeq release 57 that were used by ANNOVAR). In comparison, Ensembl provides a less curated but more comprehensive list of transcripts (n=115,901 transcripts from Ensembl version 69 that were used by ANNOVAR) based on information from the Consensus Coding Sequence (CCDS, Pruitt et al 2009), Human and Vertebrate Analysis and Annotation (HAVANA), Vertebrate Genome Annotation (Vega, Ashurst et al 2005), ENCODE (The ENCODE Project Consortium 2012) and GENCODE (Searle et al 2010).

ANNOVAR is a rapid, efficient tool to annotate functional consequences of genetic variation from high-throughput sequencing data. wANNOVAR provides easy and intuitive web-based access to the most popular functionalities of the ANNOVAR software

This tool will annotate variants using specified gene annotations, regions, and filtering databases. Input is a VCF dataset, and output is a table of annotations for each variant in the VCF dataset or a VCF dataset with the annotations in INFO fields.

ANNOVAR Website: http://annovar.openbioinformatics.org
Paper: http://nar.oxfordjournals.org/content/38/16/e164

### Download Annovar

Annovar is not opensource and we could not install it with the conda command. However, it is free for academic use.
Please register on the Annovar website and request a copy for this workshop http://download.openbioinformatics.org/annovar_download_form.php
After registering you will recive an email with a link to download the source code of Annovar on your personal computer. Please click on the link in the email, download the Annovar source code and upload on your Openstack instance using FileZilla (https://filezilla-project.org).
We have used FileZilla before in a previous part of the practical (https://github.com/KHP-Informatics/NGS_Workshop_Advanced_Bioinformatics/tree/master/assessing_quality), please see the "assessing_quality" workshop if you do not remember how to use FileZilla to download/upload files to/from your OpenStack instance.

After uploading the annovar.latest.tar.gz file onto OpenStack, please unpack it and set annovar up:

'''
$ tar -zxvf annovar.latest.tar.gz

'''

### Download Annovar db

```bash
$ ./annovar_db_download.sh
```

### VCF to Annovar input format

```bash
$ convert2annovar.pl -format vcf4 gatk.vcf > HG00403.chr20.gatk.avinput
```
```bash
NOTICE: Finished reading 3481 lines from VCF file
NOTICE: A total of 3443 locus in VCF file passed QC threshold, representing 3105 SNPs (2069 transitions and 1036 transversions) and 338 indels/substitutions
NOTICE: Finished writing 3105 SNP genotypes (2069 transitions and 1036 transversions) and 338 indels/substitutions for 1 sample
```

### Run Annovar table function

* script
```bash
$ perl annovar.pl -i HG00403.chr20.gatk.avinput -r hg38 -o HG00403.chr20.gatk
```
* tab output
```bash
$ table_annovar.pl HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -buildver hg38 -out HG00403.chr20.gatk -remove -protocol refGene,ensGene,cytoBand,genomicSuperDups,gwasCatalog,avsnp150,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur,1000g2015aug_eas,1000g2015aug_sas,nci60,cosmic85_coding,clinvar_20180603,gnomad_genome,exac03,intervar_20180118,dbnsfp31a_interpro -operation g,g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -otherinfo -nastring NA
```
* csv output
```bash
$ table_annovar.pl HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -buildver hg38 -out HG00403.chr20.gatk -remove -protocol refGene,ensGene,cytoBand,genomicSuperDups,gwasCatalog,avsnp150,esp6500siv2_all,1000g2015aug_all,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eur,1000g2015aug_eas,1000g2015aug_sas,nci60,cosmic85_coding,clinvar_20180603,gnomad_genome,exac03,intervar_20180118,dbnsfp31a_interpro, -operation g,g,r,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -otherinfo -nastring . -csvout
```

### Easy Run Annovar talbe function
```bash
$ perl annovar.pl -i demo_sample.avinput -r hg19 -o demo_sample
```

* **Gene-based (g)**
    * refGene
    * ensGene

* **Region-based (r)**
    * cytoBand
    * genomicSuperDups
    * gwasCatalog

* **Filter-based (f)**
    * avsnp150
    * esp6500siv2_all
    * 1000g2015aug_all
    * 1000g2015aug_afr
    * 1000g2015aug_amr
    * 1000g2015aug_eur
    * 1000g2015aug_eas
    * 1000g2015aug_sas
    * gnomAD genome
    * nci60
    * cosmic85
    * clinvar_20180603
    * exac03
    * intervar_20180118
    * dbnsfp31a_interpro


### Processing Information
```bash
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=refGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg38 -dbtype refGene -outfile HG00403.chr20.gatk.refGene -exonsort HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Output files were written to HG00403.chr20.gatk.refGene.variant_function, HG00403.chr20.gatk.refGene.exonic_variant_function
NOTICE: the queryfile contains 3443 lines
NOTICE: threading is disabled for gene-based annotation on file with less than 1000000 input lines
NOTICE: Reading gene annotation from /home/workshop/ngs_tools/annovar/humandb/hg38_refGene.txt ... Done with 71041 transcripts (including 17412 without coding sequence annotation) for 27813 unique genes
NOTICE: Processing next batch with 3443 unique variants in 3443 input lines
NOTICE: Reading FASTA sequences from /home/workshop/ngs_tools/annovar/humandb/hg38_refGeneMrna.fa ... Done with 7 sequences
WARNING: A total of 515 sequences will be ignored due to lack of correct ORF annotation
-----------------------------------------------------------------
NOTICE: Processing operation=g protocol=ensGene

NOTICE: Running with system command <annotate_variation.pl -geneanno -buildver hg38 -dbtype ensGene -outfile HG00403.chr20.gatk.ensGene -exonsort HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Output files were written to HG00403.chr20.gatk.ensGene.variant_function, HG00403.chr20.gatk.ensGene.exonic_variant_function
NOTICE: the queryfile contains 3443 lines
NOTICE: threading is disabled for gene-based annotation on file with less than 1000000 input lines
NOTICE: Reading gene annotation from /home/workshop/ngs_tools/annovar/humandb/hg38_ensGene.txt ... Done with 89732 transcripts (including 28806 without coding sequence annotation) for 42087 unique genes
NOTICE: Processing next batch with 3443 unique variants in 3443 input lines
NOTICE: Reading FASTA sequences from /home/workshop/ngs_tools/annovar/humandb/hg38_ensGeneMrna.fa ... Done with 13 sequences
WARNING: A total of 361 sequences will be ignored due to lack of correct ORF annotation
-----------------------------------------------------------------
NOTICE: Processing operation=r protocol=cytoBand

NOTICE: Running with system command <annotate_variation.pl -regionanno -dbtype cytoBand -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Output file is written to HG00403.chr20.gatk.hg38_cytoBand
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_cytoBand.txt ... Done with 1293 regions
NOTICE: Finished region-based annotation on 861 genetic variants
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_cytoBand.txt ... Done with 1293 regions
NOTICE: Finished region-based annotation on 861 genetic variants
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_cytoBand.txt ... Done with 1293 regions
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_cytoBand.txt ... Done with 1293 regions
NOTICE: Finished region-based annotation on 860 genetic variants
NOTICE: Finished region-based annotation on 861 genetic variants
-----------------------------------------------------------------
NOTICE: Processing operation=r protocol=genomicSuperDups

NOTICE: Running with system command <annotate_variation.pl -regionanno -dbtype genomicSuperDups -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Output file is written to HG00403.chr20.gatk.hg38_genomicSuperDups
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_genomicSuperDups.txt ... Done with 69894 regions
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_genomicSuperDups.txt ... Done with 69894 regions
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_genomicSuperDups.txt ... Done with 69894 regions
NOTICE: Finished region-based annotation on 861 genetic variants
NOTICE: Finished region-based annotation on 861 genetic variants
NOTICE: Finished region-based annotation on 861 genetic variants
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_genomicSuperDups.txt ... Done with 69894 regions
NOTICE: Finished region-based annotation on 860 genetic variants
-----------------------------------------------------------------
NOTICE: Processing operation=r protocol=gwasCatalog

NOTICE: Running with system command <annotate_variation.pl -regionanno -dbtype gwasCatalog -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Output file is written to HG00403.chr20.gatk.hg38_gwasCatalog
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_gwasCatalog.txt ... Done with 84414 regions
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_gwasCatalog.txt ... Done with 84414 regions
NOTICE: Finished region-based annotation on 861 genetic variants
NOTICE: Finished region-based annotation on 861 genetic variants
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_gwasCatalog.txt ... Done with 84414 regions
NOTICE: Reading annotation database /home/workshop/ngs_tools/annovar/humandb/hg38_gwasCatalog.txt ... Done with 84414 regions
NOTICE: Finished region-based annotation on 860 genetic variants
NOTICE: Finished region-based annotation on 861 genetic variants
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=avsnp150

NOTICE: Running system command <annotate_variation.pl -filter -dbtype avsnp150 -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_avsnp150_dropped, other variants are written to HG00403.chr20.gatk.hg38_avsnp150_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 28304406 and the number of bins to be scanned is 789
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_avsnp150.txt...Done
NOTICE: Database index loaded. Total number of bins is 28304406 and the number of bins to be scanned is 766
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_avsnp150.txt...Done
NOTICE: Database index loaded. Total number of bins is 28304406 and the number of bins to be scanned is 753
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_avsnp150.txt...Done
NOTICE: Database index loaded. Total number of bins is 28304406 and the number of bins to be scanned is 777
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_avsnp150.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=esp6500siv2_all

NOTICE: Running system command <annotate_variation.pl -filter -dbtype esp6500siv2_all -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: the --dbtype esp6500siv2_all is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_esp6500siv2_all_dropped, other variants are written to HG00403.chr20.gatk.hg38_esp6500siv2_all_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 683825 and the number of bins to be scanned is 1
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_esp6500siv2_all.txt...Done
NOTICE: Database index loaded. Total number of bins is 683825 and the number of bins to be scanned is 6
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_esp6500siv2_all.txt...Done
NOTICE: Database index loaded. Total number of bins is 683825 and the number of bins to be scanned is 19
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_esp6500siv2_all.txt...Done
NOTICE: Database index loaded. Total number of bins is 683825 and the number of bins to be scanned is 9
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_esp6500siv2_all.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=1000g2015aug_all

NOTICE: Running system command <annotate_variation.pl -filter -dbtype 1000g2015aug_all -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_ALL.sites.2015_08_dropped, other variants are written to HG00403.chr20.gatk.hg38_ALL.sites.2015_08_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 2821635 and the number of bins to be scanned is 498
NOTICE: Database index loaded. Total number of bins is 2821635 and the number of bins to be scanned is 492
NOTICE: Database index loaded. Total number of bins is 2821635 and the number of bins to be scanned is 503
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_ALL.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_ALL.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_ALL.sites.2015_08.txt...Done
NOTICE: Database index loaded. Total number of bins is 2821635 and the number of bins to be scanned is 516
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_ALL.sites.2015_08.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=1000g2015aug_afr

NOTICE: Running system command <annotate_variation.pl -filter -dbtype 1000g2015aug_afr -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_AFR.sites.2015_08_dropped, other variants are written to HG00403.chr20.gatk.hg38_AFR.sites.2015_08_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 2817370 and the number of bins to be scanned is 491
NOTICE: Database index loaded. Total number of bins is 2817370 and the number of bins to be scanned is 503
NOTICE: Database index loaded. Total number of bins is 2817370 and the number of bins to be scanned is 497
NOTICE: Database index loaded. Total number of bins is 2817370 and the number of bins to be scanned is 516
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_AFR.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_AFR.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_AFR.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_AFR.sites.2015_08.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=1000g2015aug_amr

NOTICE: Running system command <annotate_variation.pl -filter -dbtype 1000g2015aug_amr -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_AMR.sites.2015_08_dropped, other variants are written to HG00403.chr20.gatk.hg38_AMR.sites.2015_08_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 2812567 and the number of bins to be scanned is 497
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_AMR.sites.2015_08.txt...Done
NOTICE: Database index loaded. Total number of bins is 2812567 and the number of bins to be scanned is 503
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_AMR.sites.2015_08.txt...Done
NOTICE: Database index loaded. Total number of bins is 2812567 and the number of bins to be scanned is 490
NOTICE: Database index loaded. Total number of bins is 2812567 and the number of bins to be scanned is 516
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_AMR.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_AMR.sites.2015_08.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=1000g2015aug_eur

NOTICE: Running system command <annotate_variation.pl -filter -dbtype 1000g2015aug_eur -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_EUR.sites.2015_08_dropped, other variants are written to HG00403.chr20.gatk.hg38_EUR.sites.2015_08_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 2809469 and the number of bins to be scanned is 516
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_EUR.sites.2015_08.txt...Done
NOTICE: Database index loaded. Total number of bins is 2809469 and the number of bins to be scanned is 498
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_EUR.sites.2015_08.txt...Done
NOTICE: Database index loaded. Total number of bins is 2809469 and the number of bins to be scanned is 503
NOTICE: Database index loaded. Total number of bins is 2809469 and the number of bins to be scanned is 492
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_EUR.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_EUR.sites.2015_08.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=1000g2015aug_eas

NOTICE: Running system command <annotate_variation.pl -filter -dbtype 1000g2015aug_eas -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_EAS.sites.2015_08_dropped, other variants are written to HG00403.chr20.gatk.hg38_EAS.sites.2015_08_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 2810247 and the number of bins to be scanned is 516
NOTICE: Database index loaded. Total number of bins is 2810247 and the number of bins to be scanned is 503
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_EAS.sites.2015_08.txt...Done
NOTICE: Database index loaded. Total number of bins is 2810247 and the number of bins to be scanned is 497
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_EAS.sites.2015_08.txt...Done
NOTICE: Database index loaded. Total number of bins is 2810247 and the number of bins to be scanned is 491
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_EAS.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_EAS.sites.2015_08.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=1000g2015aug_sas

NOTICE: Running system command <annotate_variation.pl -filter -dbtype 1000g2015aug_sas -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_SAS.sites.2015_08_dropped, other variants are written to HG00403.chr20.gatk.hg38_SAS.sites.2015_08_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 2813352 and the number of bins to be scanned is 516
NOTICE: Database index loaded. Total number of bins is 2813352 and the number of bins to be scanned is 491
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_SAS.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_SAS.sites.2015_08.txt...Done
NOTICE: Database index loaded. Total number of bins is 2813352 and the number of bins to be scanned is 503
NOTICE: Database index loaded. Total number of bins is 2813352 and the number of bins to be scanned is 496
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_SAS.sites.2015_08.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_SAS.sites.2015_08.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=nci60

NOTICE: Running system command <annotate_variation.pl -filter -dbtype nci60 -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: the --dbtype nci60 is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_nci60_dropped, other variants are written to HG00403.chr20.gatk.hg38_nci60_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 81459 and the number of bins to be scanned is 24
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_nci60.txt...Done
NOTICE: Database index loaded. Total number of bins is 81459 and the number of bins to be scanned is 1
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_nci60.txt...Done
NOTICE: Database index loaded. Total number of bins is 81459 and the number of bins to be scanned is 14
NOTICE: Database index loaded. Total number of bins is 81459 and the number of bins to be scanned is 13
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_nci60.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_nci60.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=cosmic85_coding

NOTICE: Running system command <annotate_variation.pl -filter -dbtype cosmic85_coding -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4>
NOTICE: the --dbtype cosmic85_coding is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_cosmic85_coding_dropped, other variants are written to HG00403.chr20.gatk.hg38_cosmic85_coding_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_cosmic85_coding.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_cosmic85_coding.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_cosmic85_coding.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_cosmic85_coding.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=clinvar_20180603
NOTICE: Finished reading 5 column headers for '-dbtype clinvar_20180603'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype clinvar_20180603 -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4 -otherinfo>
NOTICE: the --dbtype clinvar_20180603 is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_clinvar_20180603_dropped, other variants are written to HG00403.chr20.gatk.hg38_clinvar_20180603_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Database index loaded. Total number of bins is 44007 and the number of bins to be scanned is 12
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_clinvar_20180603.txt...Done
NOTICE: Database index loaded. Total number of bins is 44007 and the number of bins to be scanned is 2
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_clinvar_20180603.txt...Done
NOTICE: Database index loaded. Total number of bins is 44007 and the number of bins to be scanned is 15
NOTICE: Database index loaded. Total number of bins is 44007 and the number of bins to be scanned is 26
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_clinvar_20180603.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_clinvar_20180603.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=gnomad_genome
NOTICE: Finished reading 8 column headers for '-dbtype gnomad_genome'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype gnomad_genome -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4 -otherinfo>
NOTICE: the --dbtype gnomad_genome is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_gnomad_genome_dropped, other variants are written to HG00403.chr20.gatk.hg38_gnomad_genome_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 28084439 and the number of bins to be scanned is 753
NOTICE: Database index loaded. Total number of bins is 28084439 and the number of bins to be scanned is 789
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_gnomad_genome.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_gnomad_genome.txt...Done
NOTICE: Database index loaded. Total number of bins is 28084439 and the number of bins to be scanned is 776
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_gnomad_genome.txt...Done
NOTICE: Database index loaded. Total number of bins is 28084439 and the number of bins to be scanned is 766
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_gnomad_genome.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=exac03
NOTICE: Finished reading 8 column headers for '-dbtype exac03'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype exac03 -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4 -otherinfo>
NOTICE: the --dbtype exac03 is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_exac03_dropped, other variants are written to HG00403.chr20.gatk.hg38_exac03_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 749044 and the number of bins to be scanned is 0
NOTICE: Database index loaded. Total number of bins is 749044 and the number of bins to be scanned is 7
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_exac03.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_exac03.txt...Done
NOTICE: Database index loaded. Total number of bins is 749044 and the number of bins to be scanned is 8
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_exac03.txt...Done
NOTICE: Database index loaded. Total number of bins is 749044 and the number of bins to be scanned is 15
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_exac03.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=intervar_20180118
NOTICE: Finished reading 29 column headers for '-dbtype intervar_20180118'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype intervar_20180118 -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4 -otherinfo>
NOTICE: the --dbtype intervar_20180118 is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_intervar_20180118_dropped, other variants are written to HG00403.chr20.gatk.hg38_intervar_20180118_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 535845 and the number of bins to be scanned is 8
NOTICE: Database index loaded. Total number of bins is 535845 and the number of bins to be scanned is 0
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_intervar_20180118.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_intervar_20180118.txt...Done
NOTICE: Database index loaded. Total number of bins is 535845 and the number of bins to be scanned is 6
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_intervar_20180118.txt...Done
NOTICE: Database index loaded. Total number of bins is 535845 and the number of bins to be scanned is 3
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_intervar_20180118.txt...Done
-----------------------------------------------------------------
NOTICE: Processing operation=f protocol=dbnsfp31a_interpro
NOTICE: Finished reading 1 column headers for '-dbtype dbnsfp31a_interpro'

NOTICE: Running system command <annotate_variation.pl -filter -dbtype dbnsfp31a_interpro -buildver hg38 -outfile HG00403.chr20.gatk HG00403.chr20.gatk.avinput /home/workshop/ngs_tools/annovar/humandb/ -thread 4 -otherinfo>
NOTICE: the --dbtype dbnsfp31a_interpro is assumed to be in generic ANNOVAR database format
NOTICE: Variants matching filtering criteria are written to HG00403.chr20.gatk.hg38_dbnsfp31a_interpro_dropped, other variants are written to HG00403.chr20.gatk.hg38_dbnsfp31a_interpro_filtered
NOTICE: the queryfile contains 3443 lines
NOTICE: Creating new threads for query line 1 to 861
NOTICE: Creating new threads for query line 862 to 1722
NOTICE: Creating new threads for query line 1723 to 2583
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Creating new threads for query line 2584 to 3443
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 861 unique variants in 861 input lines
NOTICE: Processing next batch with 860 unique variants in 860 input lines
NOTICE: Database index loaded. Total number of bins is 275396 and the number of bins to be scanned is 0
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_dbnsfp31a_interpro.txt...Done
NOTICE: Database index loaded. Total number of bins is 275396 and the number of bins to be scanned is 2
NOTICE: Database index loaded. Total number of bins is 275396 and the number of bins to be scanned is 2
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_dbnsfp31a_interpro.txt...Done
NOTICE: Database index loaded. Total number of bins is 275396 and the number of bins to be scanned is 4
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_dbnsfp31a_interpro.txt...Done
NOTICE: Scanning filter database /home/workshop/ngs_tools/annovar/humandb/hg38_dbnsfp31a_interpro.txt...Done
-----------------------------------------------------------------
NOTICE: Multianno output file is written to HG00403.chr20.gatk.hg38_multianno.txt
Finished
```

---

### Database download log
```bash
./annovar_db_download.sh
```

```bash
Create humandb folder done!
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_refGene.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_refGeneMrna.fa.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_refGeneVersion.txt.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_refGene.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_refGeneMrna.fa.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_refGeneVersion.txt.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_knownGene.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_kgXref.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_knownGeneMrna.fa.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_knownGene.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_kgXref.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_knownGeneMrna.fa.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_ensGene.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_ensGeneMrna.fa.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_ensGene.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_ensGeneMrna.fa.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
--2018-08-06 17:58:27--  http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/cytoBand.txt.gz
Resolving hgdownload.cse.ucsc.edu (hgdownload.cse.ucsc.edu)... 128.114.119.163
Connecting to hgdownload.cse.ucsc.edu (hgdownload.cse.ucsc.edu)|128.114.119.163|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 9882 (9.7K) [application/x-gzip]
Saving to: ‘cytoBand.txt.gz’

cytoBand.txt.gz                              100%[============================================================================================>]   9.65K  --.-KB/s    in 0s

2018-08-06 17:58:27 (204 MB/s) - ‘cytoBand.txt.gz’ saved [9882/9882]

NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
--2018-08-06 17:58:34--  http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/genomicSuperDups.txt.gz
Resolving hgdownload.cse.ucsc.edu (hgdownload.cse.ucsc.edu)... 128.114.119.163
Connecting to hgdownload.cse.ucsc.edu (hgdownload.cse.ucsc.edu)|128.114.119.163|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 4202400 (4.0M) [application/x-gzip]
Saving to: ‘genomicSuperDups.txt.gz’

genomicSuperDups.txt.gz                      100%[============================================================================================>]   4.01M   685KB/s    in 7.0s

2018-08-06 17:58:41 (583 KB/s) - ‘genomicSuperDups.txt.gz’ saved [4202400/4202400]

NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/phastConsElements100way.txt.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
--2018-08-06 17:59:04--  http://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/phastConsElements100way.txt.gz
Resolving hgdownload.cse.ucsc.edu (hgdownload.cse.ucsc.edu)... 128.114.119.163
Connecting to hgdownload.cse.ucsc.edu (hgdownload.cse.ucsc.edu)|128.114.119.163|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 90172140 (86M) [application/x-gzip]
Saving to: ‘phastConsElements100way.txt.gz’

phastConsElements100way.txt.gz               100%[============================================================================================>]  85.99M  14.4MB/s    in 18s

2018-08-06 17:59:22 (4.86 MB/s) - ‘phastConsElements100way.txt.gz’ saved [90172140/90172140]

NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_1000g2015aug.zip ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_1000g2015aug.zip ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_avsnp150.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_avsnp150.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_avsnp150.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_avsnp150.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_nci60.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_nci60.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_nci60.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_nci60.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_clinvar_20180603.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_clinvar_20180603.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_clinvar_20180603.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_clinvar_20180603.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_gnomad_genome.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_gnomad_genome.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_gnomad_genome.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_gnomad_genome.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_gnomad_exome.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_gnomad_exome.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_gnomad_exome.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_gnomad_exome.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_exac03.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_exac03.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_exac03.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_exac03.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_esp6500siv2_all.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_esp6500siv2_all.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_esp6500siv2_all.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_esp6500siv2_all.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_intervar_20180118.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_intervar_20180118.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_intervar_20180118.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_intervar_20180118.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the '/home/philippe/humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_dbnsfp31a_interpro.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg19_dbnsfp31a_interpro.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the 'humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp31a_interpro.txt.gz ... OK
NOTICE: Downloading annotation database http://www.openbioinformatics.org/annovar/download/hg38_dbnsfp31a_interpro.txt.idx.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the 'humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/gwasCatalog.txt.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg19 build version, with files saved at the 'humandb' directory
NOTICE: Web-based checking to see whether ANNOVAR new version is available ... Done
NOTICE: Downloading annotation database http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/gwasCatalog.txt.gz ... OK
NOTICE: Uncompressing downloaded files
NOTICE: Finished downloading annotation files for hg38 build version, with files saved at the 'humandb' directory
```
---

* Reference by:
    * [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/)
    * [Using VAAST to Identify an X-Linked Disorder Resulting in Lethality in Male Infants Due to N-Terminal Acetyltransferase Deficiency](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3135802/)
    * [ANNOVAR: functional annotation of genetic variants from high-throughput sequencing data](http://www.ncbi.nlm.nih.gov/pubmed/20601685)
    * [Genomic variant annotation and prioritization with ANNOVAR and wANNOVAR](http://www.nature.com/nprot/journal/v10/n10/full/nprot.2015.105.html)
    * [NAA10 / rs387906701](http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=rs387906701)
    * [SNPedia](http://snpedia.com/index.php/Rs387906701)
    * [NAA10 mutation causing a novel intellectual disability syndrome with Long QT due to N-terminal acetyltransferase impairment](http://www.nature.com/articles/srep16022)
    * [OMIM - Ogden Syndrome](http://www.omim.org/entry/300855)


* All Information 2018 Genomic Epidemiology Workshop Use
* Edit by [Philippe](http://github.com/geniusphil)
