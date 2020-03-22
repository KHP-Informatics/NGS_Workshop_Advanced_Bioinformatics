---
title: "Automated bash pipeline"
author: "Alfredo Iacoangeli"
date: "16/03/2020"
---

IMPORTANT: the paths to files and directories in this workshop are only examples. You will need use your own so please try not to copy and paste the commands but write them yourself matching the correct locations on your computer

## Learning Objectives:
* Learn how to automate your pipeline using a bash script 
* Use a `for loop` to run the pipeline on multiple files


Suppose you want to run a few different commands in succession – this is called “automating tasks”. See http://xkcd.com/1205/

Key caveat: none of the commands require any user live input (‘yes’, ‘no’, etc)

### New project directory tree

First let's create a directory tree for this project:

Make a directory for the DNA-seq automated pipeline analysis within the `ngs_course` folder using the `mkdir` command

```
$ cd ~/ngs_course

$ mkdir dnaseq_pipeline
```

Next you want to set up the following structure within your project directory to keep files organized:

```
dnaseq_pipeline/
  ├── data/
  ├── meta/
  ├── results/
  ├── logs/
  └── scripts/

```

Let's create a directory for our project by changing into `dnaseq` and then using `mkdir` to create the four directories.

```
$ cd ngs_course/dnaseq_pipeline

$ mkdir data meta results logs scripts
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
drwxrwsr-x 2 rsk27 rsk27   0 Jun 17 11:21 scripts/
```
Now we will create a subdirectory to to store the data we want to process. We will start with the trimmed data from the previous project "dnaseq": 

```
$ cd ~/ngs_course/dnaseq_pipeline/data

$ mkdir trimmed_fastq
```
    
The fastq data we will be working with need to be copied from dnaseq first (this include the sequencing data as well as the bed file and the reference genome:

#### IMPORTANT: please keep in mind that the total available storage on your OpenStack instance is 40 Gigabytes. You can check the total amount of storage you have used by running the command `sudo du -sh /`. If you realize that you do not have enough space to copy the data in the new project, consider moving the data instead using the `mv` command instead of `cp`

```
$ cp ~/ngs_course/dnaseq/data/trimmed_fastq/TRIMMED_DATA_FILENAME_R1 ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq

$ cp ~/ngs_course/dnaseq/data/trimmed_fastq/TRIMMED_DATA_FILENAME_R2 ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq

$ cp ~/ngs_course/dnaseq/data/chr22.genes.hg19.bed ~/ngs_course/dnaseq_pipeline/data/

$ mkdir ~/ngs_course/dnaseq/data/dnaseq_pipeline/reference

$ cp ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq/data/dnaseq_pipeline/reference/
```

### Create an example bash script

* Figure out what commands you want to run.
* Put them in a file using a text editor, like nano or pico, or TextEdit, or TextWrangler, or vi, or emacs.
* Run the file.
* For example, try putting the following shell commands in a bash script called ‘test.sh’:

```
$ for dirname in $(ls ~/ngs_course/dnaseq/)
> do
>   echo The following is the content of $dirname :
>   ls ~/ngs_course/dnaseq/$dirname
>   echo
> done
```
This commands list all the files contained in the dnaseq folders.

Now create test.sh and add the commands to the file then run it and see the poutput:

```
$ cd ~/ngs_course/dnaseq_pipeline/scripts

$ vim test.sh
```
Change to insert mode, then type in the following lines in the test.sh file:

```
#!/bin/bash #

for dirname in $(ls ~/ngs_course/dnaseq/)
do
echo The following is the content of $dirname
ls ~/ngs_course/dnaseq/$dirname
echo
done
```

Save the file and exit vim. Now let's run the new script we have created. To run a shell script you usually use the bash or sh command.



```
$ bash test.sh

The following is the content of data
chr22.genes.hg19.bed  other  reference	trimmed_fastq  untrimmed_fastq

The following is the content of logs
fastqc_summaries.txt

The following is the content of meta

The following is the content of results
fastqc_untrimmed_reads something_else 
```

