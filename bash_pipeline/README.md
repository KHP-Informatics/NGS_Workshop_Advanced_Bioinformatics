---
title: "Automated bash pipeline"
author: "Alfredo Iacoangeli"
date: "16/03/2020"
---

IMPORTANT: the paths to files and directories in this workshop are only examples. You will need use your own so please try not to copy and paste the commands but write them yourself matching the correct locations on your computer

## Learning Objectives:
* Learn how to automate your pipeline using a bash script 


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

$ mkdir untrimmed_fastq
```
    
The fastq data we will be working with need to be copied from dnaseq first (this include the sequencing data as well as the bed file and the reference genome:

#### IMPORTANT: please keep in mind that the total available storage on your OpenStack instance is 40 Gigabytes. You can check the total amount of storage you have used by running the command `sudo du -sh /`. If you realize that you do not have enough space to copy the data in the new project, consider moving the data instead using the `mv` command instead of `cp`

```
$ cp ~/ngs_course/dnaseq/data/untrimmed_fastq/WES01_chr22m_R1.fastq.gz ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq

$ cp ~/ngs_course/dnaseq/data/trimmed_fastq/WES01_chr22m_R2.fastq.gz ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq

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
>   echo                                             # we use this command to add a newline
> done
```
This commands list all the files contained in the dnaseq folders.

Now create test.sh and add the commands to the file then run it and see the output:

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
## Bash variables (repetition from the loops_scripts workshop) 
A *variable* is a common concept shared by many programming languages. Variables are essentially a symbolic/temporary name for, or a reference to, some information. Variables are analogous to "buckets", where information can be stored, maintained and modified without too much hassle. 

Extending the bucket analogy: the bucket has a name associated with it, i.e. the name of the variable, and when referring to the information in the bucket, we use the name of the bucket, and do not directly refer to the actual data stored in it.

Let's start with a simple variable that has a single number stored in it:

	$ num=25

*How do we know that we actually created the bash variable?* We can use the `echo` command to print to terminal:

	$ echo num


What do you see in the terminal? The `echo` utility takes what arguments you provide and prints to terminal. In this case it interpreted `num` as a character string and simply printed it back to us. This is because **when using the variable as an argument to the `echo` command, we explicitly use a `$` in front of it**:

	$ echo $num

Now you should see the number 25 returned to you. Did you notice that when we created the variable we just typed in the variable name? This is standard shell notation (syntax) for defining and using variables. When defining the variable (i.e. setting the value) you can just type it as is, but when **retrieving the value of a variable don't forget the `$`!** 

Variables can also store a string of character values. In the example below, we define a variable or a 'bucket' called `file`. We will put a filename `WES01_chr22m_R1.fastq.gz` as the value inside the bucket.

	$ file=WES01_chr22m_R1.fastq.gz

Once you press return, you should be back at the command prompt. Let's check what's stored inside `file`:

	$ echo $file


Let's try another command using the variable that we have created. In the last lesson, we introduced the `wc -l` command which allows us to count the number of lines in a file. We can count the number of lines in `WES01_chr22m_R1.fastq.gz` by referencing the `file` variable, but first move into the `raw_fastq` directory:

	$ cd ~/ngs_course/dnaseq/data/untrimmed_fastq
	$ wc -l $file

> *NOTE:* The variables we create in a session are system-wide, and independent of where you are in the filesystem. This is why we can reference it from any directory. However, it is only available for your current session. If you exit the cluster and login again at a later time, the variables you have created will no longer exist.

Ok, so we know variables are like buckets, and so far we have seen that bucket filled with a single value. **Variables can store more than just a single value.** They can store multiple values and in this way can be useful to carry out many things at once. Let's create a new variable called `filenames` and this time we will store *all of the filenames* in the `untrimmed_fastq` directory as values. 

To list all the filenames in the directory that have a `.fq` extension, we know the command is:

	$ ls *.fastq.gz
	
Now we want to *assign* the output of `ls` to the variable. We will give that variable the name `filenames`:

	$ filenames=`ls *.fastq.gz`

Check and see what's stored inside our newly created variable using `echo`:
	
	$ echo $filenames

Let's try the `wc -l` command again, but this time using our new variable `filenames` as the argument:

	$ wc -l $filenames
	
What just happened? Because our variable contains multiple values, the shell runs the command on each value stored in `filenames` and prints the results to screen. 

### Special Variable Types: Positional parameters

In any bash script the variables $0, $1, $2, $3, etc are the arguments passed to the script from the command line. $0 is the name of the script itself, $1 is the first argument, $2 the second, $3 the third, and so forth. [2] After $9, the arguments must be enclosed in brackets, for example, ${10}, ${11}, ${12}. See following example:

Please open a new file with vim, write the following lines into it, save it and exit:

```
$ cd ~/ngs_course/dnaseq_pipeline/scripts
$ vim test_1.sh
```
Write the following lines into:

```
#!/bin/bash

echo the name of the script is $0
echo 
echo the first argument was $1 and the second $2

```

Now run the script with your own arguments:

```
$ bash test_1.sh something "something else" 
```

Please note that we used the `" "` for the second argument because without them bash would have read "something" as the second argument and else as the third.

## A script to run Trimmomatic on the sequencing data

Now that we know  how to write simple bash scripts and run Trimmomatic on PE seq data we can write a script that takes as argument the sequencing data file names and performs data trimming. We are going to create the file pipeline.sh and add the commands lines to perform data trimming as we did in the ["trimming"](https://github.com/KHP-Informatics/NGS_Workshop_Advanced_Bioinformatics/tree/master/trimming) workshop:

```
$ vim pipeline.sh
```

Let's add the following lines to the script (please note that you need to modify the command line according to you directory tree):

```
#!/bin/bash

trimmomatic PE  \
  -threads 4 \
  -phred33 \
  $1 $2 \
  -baseout ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_fastq/trimmed_data \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50
```
Now let's run the script:

```
$ bash pipeline.sh ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R1.fastq.gz \ 				~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R2.fastq.gz
```
Any errors or standard output from the commands run within the script will be directed onto the terminal. If not error is produced and the terminal shows something like the following, please go to your "trimmed_fastq" directory and see that the trimmed data was effectively generated.

```
trimmomaticPE: Started with arguments:
 -threads 4 -phred33 WES01_chr22m_R1.fastq WES01_chr22m_R2.fastq -baseout /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_fastq/trimmed_data ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50
Using templated Output files: /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_fastq/trimmed_data_1P /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_fastq/trimmed_data_1U /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_fastq/trimmed_data_2P /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_fastq/trimmed_data_2U
Using PrefixPair: 'AGATGTGTATAAGAGACAG' and 'AGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTCCGAGCCCACGAGAC'
Using Long Clipping Sequence: 'CTGTCTCTTATACACATCTGACGCTGCCGACGA'
ILLUMINACLIP: Using 1 prefix pairs, 4 forward/reverse sequences, 0 forward only sequences, 0 reverse only sequences
Input Read Pairs: 505987 Both Surviving: 432839 (85.54%) Forward Only Surviving: 67687 (13.38%) Reverse Only Surviving: 2163 (0.43%) Dropped: 3298 (0.65%)
TrimmomaticPE: Completed successfully

$ cd /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_fastq/

$ ls

trimmed_data_1P  trimmed_data_1U  trimmed_data_2P  trimmed_data_2U
```

## Add the FastQC analysis to your script

Now that we have generated a pipeline.sh script that can perform data trimming on any input fastq files, we want the script to run FastQC on the generated trimmed data. Open pipeline.sh with vim and add the following command line:

```

fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
	/home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P
	
mkdir ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads

mv ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/*fastqc* ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads/

```

Now run pipeline.sh and see if it generates the FastQC reports in the "fastqc_trimmed_reads" directory. Please note that every time you run pipeline.sh the script will run all its commands and this can be quite time-consuming and risky in case it overwrites something important. For the aim of this practical exercise this is not important but please keep this in mind when you will make another automated pipeline in the future.

```
$ bash pipeline.sh ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R1.fastq.gz \ 				~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R2.fastq.gz

$ cd ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads/

$ ls
trimmed_data_1P_fastqc.html  trimmed_data_1P_fastqc.zip  trimmed_data_2P_fastqc.html  trimmed_data_2P_fastqc.zip
```

## Final excercise

Now that you know how to generate a bash script to run a sequence of commands, please complete the pipeline.sh with all the steps of a standard DNA-seq analysis pipeline by adding the commands lines from the other workshops to pipeline.sh and run the whole analysis pipeline. The script should be able to perform the following steps automatically by running only one command:

```
bash pipeline.sh ~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R1.fastq.gz \ 				~/ngs_course/dnaseq_pipeline/data/untrimmed_fastq/WES01_chr22m_R2.fastq.gz
```

- Trimming
- QC (fastQC)
- Alignement
- Duplicate marking
- Post-alignemtn read filtering
- Variant Calling
- Variant filtering
- Variant annotation.

For inspiration please see [the last slides of the module introduction lecture](https://github.com/KHP-Informatics/NGS_Workshop_Advanced_Bioinformatics/blob/master/bash_pipeline/exsample_bash_script.pdf) (please note that Picard and Trimmomatic command lines in the slides are different from what you are using, do not copy and paste)

## IMPORTANT: 
### 1) Please not that the pipeline.sh does not need to install the tools or download the reference or databases in this occasion. We assume that you already set up an appropriate environment in the project set up phase.
### 2) Some of the steps of the pipeline might take a little while to run, e.g. reference indexing. Please take this into account when you test it and consider creating the index beforehand as a project set up step and not as a step in the pipeline.
### 3) Please comment your script extensively by describing what each command does (you can add comments to a script by [putting the # before the comment](https://www.tutorialkart.com/bash-shell-scripting/bash-comments/) to make it clear to understand
### 4) Check that the results are the same as the ones you generated in the dnaseq project
### 5) Upload your final script to your Github repository, the one you have generated on the first day [Link](https://github.com/KHP-Informatics/NGS_Workshop_Advanced_Bioinformatics/blob/master/basic-git-bash-task.pdf)





