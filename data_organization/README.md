---
title: "How to set up a bioinformatics project"
author: "Alfredo Iacoangeli"
date: "16/03/2020"
---

Approximate time: 20 minutes

## Learning Objectives:

* Setting up your project space for an NGS workflow

## DNA-seq Workflow

The workflow we will be using for our DNA-Seq analysis is provided below with a brief description of each step.

<img src=download.png width=500>

1. Sequence quality reports - Assessing quality using FastQC
2. Sequence quality filters - Trimming and/or filtering reads (if necessary)
3. Alignment - Index the reference genome for use by BWA
4. Alignment - Align reads to reference genome using BWA 
5. Post alignment Read Filters - Duplicate marking
6. Variant calling - Variant calling using Freebayes
7. Variant quality filters - Variant filtering
8. Variant annotation - Variant annotation using Annovar

Before we jump into any kind of bioinformatics analysis, we are going to take a moment and spend some time discussing ways in which you can get organized.

## Project organization

Project organization is one of the most important parts of a sequencing project, but is **often overlooked in the excitement to get a first look at new data**. While it's best to get yourself organized before you begin analysis, it's never too late to start.

In the most important ways, the methods and approaches we need in bioinformatics are the same ones we need at the bench or in the field - **planning, documenting, and organizing** will be the key to good reproducible science. 

### Planning 

You should approach your sequencing project in a very similar way to how you do a biological experiment, and ideally, begins with **experimental design**. We're going to assume that you've already designed a beautiful sequencing experiment to address your biological question, collected appropriate samples, and that you have enough statistical power. 

### Organizing

Every computational analysis you do is going to spawn many files, and inevitability, you'll want to run some of those analysis again. Genomics projects can quickly accumulate hundreds of files across tens of folders. Before you start any analysis it is best to first get organized and **create a planned storage space for your workflow**.

We will start by creating a directory that we can use for the rest of the DNA-seq session:

First, make sure that you are in your home directory,

```
$ pwd
```
this should give the result: `/home/user_name`

**Tip** If you were not in your home directory, the easiest way to get there is to enter the command `cd` - which always returns you to home. 

Now, make a directory for the DNA-seq analysis within the `ngs_course` folder using the `mkdir` command

```
$ mkdir ngs_course

$ mkdir ngs_course/dnaseq
```

Next you want to set up the following structure within your project directory to keep files organized:

```
dnaseq/
  ├── data/
  ├── meta/
  ├── results/
  └── logs/

```

*This is a generic structure and can be tweaked based on personal preferences.* A brief description of what might be contained within the different sub-directories is provided below:

* **`data/`**: This folder is usually reserved for any raw data files that you start with. 

* **`meta/`**: This folder contains any information that describes the samples you are using, which we often refer to as metadata. 

* **`results/`**: This folder will contain the output from the different tools you implement in your workflow. To stay organized, you should create sub-folders specific to each tool/step of the workflow. 

* **`logs/`**: It is important to keep track of the commands you run and the specific pararmeters you used, but also to have a record of any standard output that is generated while running the command. 


Let's create a directory for our project by changing into `dnaseq` and then using `mkdir` to create the four directories.

```
$ cd ngs_course/dnaseq
$ mkdir data meta results logs
``` 

Verify that you have created the directories:

```
$ ls -lF
```
 
If you have created these directories, you should get the following output from that command:

```
drwxrwsr-x 2 rsk27 rsk27   0 Jun 17 11:21 data/
drwxrwsr-x 2 rsk27 rsk27   0 Jun 17 11:21 logs/
drwxrwsr-x 2 rsk27 rsk27   0 Jun 17 11:21 meta/
drwxrwsr-x 2 rsk27 rsk27   0 Jun 17 11:21 results/
```
Now we will create the subdirectories to setup for our DNA-Seq analysis, and populate them with data where we can. The first step will be checking the quality of our data, and trimming the files if necessary. We need to create two directories within the `data` directory, one folder for untrimmed reads and another for our trimmed reads: 

```
$ cd ~/ngs_course/dnaseq/data
$ mkdir untrimmed_fastq
$ mkdir trimmed_fastq
```
    
The raw_fastq data we will be working with need to be downloaded first:

```
$ wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/WES01_chr22m_R1.fastq.gz

$ wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/WES01_chr22m_R2.fastq.gz

$ wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/chr22.genes.hg19.bed

```

Then we need to copy the raw fastq files to our `untrimmed_fastq` directory and the bed file in the `data` directory:

`$ mv *fastq.gz ~/ngs_course/dnaseq/data/untrimmed_fastq`

`$ mv chr22.genes.hg19.bed ~/ngs_course/dnaseq/data`

Later in the workflow when we perform alignment, we will require reference files to map against. Please download it and then  you can copy it into `data`:

`$ wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz`

`$ mv hg19.fa.gz ~/ngs_course/dnaseq/data/`

### Documenting

For all of those steps, collecting specimens, extracting DNA, prepping your samples, you've likely kept a lab notebook that details how and why you did each step, but **documentation doesn't stop at the sequencer**! 

 
#### README

Keeping notes on what happened in what order, and what was done, is essential for reproducible research. If you don’t keep good notes, then you will forget what you did pretty quickly, and if you don’t know what you did, noone else has a chance. After setting up the filesystem and running a workflow it is useful to have a **README file within your project** directory. This file will usually contain a quick one line summary about the project and any other lines that follow will describe the files/directories found within it. Within each sub-directory you can also include README files 
to describe the analysis and the files that were generated. 


#### Log files

In your lab notebook, you likely keep track of the different reagents and kits used for a specific protocol. Similarly, recording information about the tools and and parameters is imporatant for documenting your computational experiments. 

* Keep track of software versions used
* Record information on parameters used and summary statistics at every step (e.g., how many adapters were removed, how many reads did not align)

> Different tools have different ways of reporting log messages and you might have to experiment a bit to figure out what output to capture. You can redirect standard output with the `>` symbol which is equivalent to `1> (standard out)`; other tools might require you to use `2>` to re-direct the `standard error` instead. 
 
***

**Exercise**

1. Take a moment to create a README for the `dnaseq` folder (hint: use `vim`, or `nano` or any text editor of your choice to create the file). Give a short description of the project and brief descriptions of the types of file you would be storing within each of the sub-directories. 

***

###Install the tools you need for the following workshops

Before going ahead with the other NGS practical materials we need to install the tools requiered for the analyses. In the Kotacoda workshops you have learned how to deploy Anaconda locally and use it to install other tools (https://katacoda.com/sjnewhouse/scenarios/bioconda_genmed)

Please install Anaconda locally on your OpenStack instance and use it to install Trimmomatic, Fastqc, Samtools, Picard, Bedtools, BWA, Freebayes and any other tools you will need for this module.

'''

$ cd ~/

$ wget https://repo.anaconda.com/archive/Anaconda3-5.3.0-Linux-x86_64.sh

$ chmod +x ./Anaconda3-5.3.0-Linux-x86_64.sh

$ bash ./Anaconda3-5.3.0-Linux-x86_64.sh

$ source ~/.bashrc

$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge

$ conda install samtools
$ conda install bwa
$ conda install freebayes
$ conda install picard
$ conda install bedtools
$ conda install trimmomatic
$ conda install fastqc
$ conda install any_other_tool_you_need

'''

----

*This lesson has been developed using material produced by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
