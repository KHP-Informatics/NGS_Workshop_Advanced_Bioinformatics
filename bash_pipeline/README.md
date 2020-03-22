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
    
The fastq data we will be working with need to be copied from dnaseq first:

#### IMPORTANT: please keep in mind that the total available storage on your OpenStack instance is 40 Gigabytes. You can check the total amount of storage you have used by running the command `du -sh /`. If you realize that you do not have enough space to copy the data in the new project, consider moving the data instead using the `mv` command instead of `cp`

```
$ cp ~/ngs_course/dnaseq/data/trimmed_fastq/TRIMMED_DATA_FILENAME_R1 ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq

$ cp ~/ngs_course/dnaseq/data/trimmed_fastq/TRIMMED_DATA_FILENAME_R2 ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq
```

* Figure out what commands you want to run.
* Put them in a file using a text editor, like nano or pico, or TextEdit, or TextWrangler, or vi, or emacs.
* Run the file.
* For example, try putting the following shell commands in a bash script called ‘test.sh’:
