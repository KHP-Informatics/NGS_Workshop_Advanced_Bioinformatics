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
>   echo                                             # we use this command to add a newline
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
## Bash variables
A *variable* is a common concept shared by many programming languages. Variables are essentially a symbolic/temporary name for, or a reference to, some information. Variables are analogous to "buckets", where information can be stored, maintained and modified without too much hassle. 

Extending the bucket analogy: the bucket has a name associated with it, i.e. the name of the variable, and when referring to the information in the bucket, we use the name of the bucket, and do not directly refer to the actual data stored in it.

Let's start with a simple variable that has a single number stored in it:

	$ num=25

*How do we know that we actually created the bash variable?* We can use the `echo` command to print to terminal:

	$ echo num


What do you see in the terminal? The `echo` utility takes what arguments you provide and prints to terminal. In this case it interpreted `num` as a a character string and simply printed it back to us. This is because **when using the variable as an argument to the `echo` command, we explicitly use a `$` in front of it**:

	$ echo $num

Now you should see the number 25 returned to you. Did you notice that when we created the variable we just typed in the variable name? This is standard shell notation (syntax) for defining and using variables. When defining the variable (i.e. setting the value) you can just type it as is, but when **retrieving the value of a variable don't forget the `$`!** 

Variables can also store a string of character values. In the example below, we define a variable or a 'bucket' called `file`. We will put a filename `WES01_chr22m_R1.fastq.gz` as the value inside the bucket.

	$ file=WES01_chr22m_R1.fastq.gz

Once you press return, you should be back at the command prompt. Let's check what's stored inside `file`:

	$ echo $file


Let's try another command using the variable that we have created. In the last lesson, we introduced the `wc -l` command which allows us to count the number of lines in a file. We can count the number of lines in `Mov10_oe_1.subset.fq` by referencing the `file` variable, but first move into the `raw_fastq` directory:

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
