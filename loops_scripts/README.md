---
title: "The Shell: Loops & Scripts"
author: "Alfredo Iacoangeli"
date: "16/03/2020"
---

## IMPORTANT: the paths to files and directories in this workshop are only examples. You will need use your own so please try not to copy and paste the commands but write them yourself matching the correct locations on your computer

## Learning Objectives

* Learn how to operate on multiple files 
* Capture previous commands into a script to re-run later
* Automating a workflow with scripts

Now that you've been using quite a number of commands to interrogate your data, 
wouldn't it be great if you could do this for each set of data that comes in, without having to manually re-type the commands?

Welcome to the beauty and purpose of shell scripts.

## Shell scripts

Shell scripts are text files that contain commands we want to run. As with any file, you can give a shell script any name and usually have the extension `.sh`. For historical reasons, a bunch of commands saved in a file is usually called a shell script, but make no mistake, this is actually a small program. 


We are going to take the commands we repeat frequently and save them into a file so that we can **re-run all those operations** again later by typing **one single command**. Let's write a shell script that will do two things:

1. Tell us what our current working directory is
2. Lists the contents of the directory 

First let's move into the `dnaseq` directory and open a new file using `vim`:

	$ cd ~/ngs_course/dnaseq
	$ vim listing.sh
	
Change to insert mode, then type in the following lines in the `listing.sh` file:

	echo "Your current working directory is:"
	pwd
	echo "These are the contents of this directory:"
	ls -l 

>The `echo` command is a utility for writing to standard output. By providing text in quotations after the command we indicated what it is we wanted written

Save the file and exit `vim`. Now let's run the new script we have created. To run a shell script you usually use the `bash` or `sh` command.

	$ sh listing.sh
	
> Did it work like you expected?
> 
> Were the `echo` commands helpful in letting you know what came next?

This is a very simple shell script. In this session and in upcoming sessions, we will be learning how to write more complex ones. You will see how the power of scripts can make our lives much easier.

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

> Try using some of the other commands we learned in previous lessons (i.e `head`, `tail`) on the `filename` variable. 

## Loops

Another powerful concept in the Unix shell is the concept of "Loops". We have just shown you that you can run a single command on multiple files by creating a variable whose values are the filenames that you wish to work on. But what if you want to **run a sequence of multiple commands, on multiple files**? This is where loop come in handy!

Looping is a concept shared by several programming languages, and its implementation in bash is very similar to other languages. 

The structure or the syntax of (*for*) loops in bash is as follows:

```bash
$ for (variable_name) in (list)
> do
>   (command $variable_name) 
> done
```

where the ***variable_name*** defines (or initializes) a variable that takes the value of every member of the specified ***list*** one at a time. At each iteration, the loop retrieves the value stored in the variable (which is a member of the input list) and runs through the commands indicated between the `do` and `done` one at a time. *This syntax/structure is virtually set in stone.* 

For example, we can run the same commands (`echo` and `wc -l`) used in the "Bash variables" section but this time run them sequentially on each file:

```
$ for filename in *.fastq.gz
> do
>   echo $filename
>   wc -l $filename
> done
````

#### What does this loop do? 
Most simply, it writes to the terminal (`echo`) the name of the file and the number of lines (`wc -l`) for each files that end in `.fastq.gz` in the current directory. The output is almost identical to what we had before.

In this case the list of files is specified using the asterisk wildcard: `*.fastq.gz`, i.e. all files that end in `.fastq.gz`. Then, we execute 2 commands between the `do` and `done`. With a loop, we execute these commands for each file one at a time. For each iteration the filename gets stored in the temporary variable called `filename`. Once the commands are executed for one file, the loop then stores the next filename in `filename` and executes the same commands on the next file. In the long run, it's best to use a name that will help point out a variable's function, so your future self will understand what you are thinking now.
 
Essentially, **the number of loops == the number of items in the list**, in our case that is 6 times since we have 6 files in `~/ngs_course/dnaseq/data/untrimmed_fastq/` that end in `.fastq.gz`. This is done by changing the value of the `filename` variable 6 times. 
Pretty simple and cool, huh?

## Automating with Scripts
	
Now that you've learned how to use loops and variables, let's put this processing power to work. Imagine, if you will, a series of commands that would do the following for us each time we get a new data set:

- Use for loop to iterate over each FASTQ file
- Dump out bad reads into a new file
- Get the count of the number of bad reads
- And write the filename and number of bad reads to a summary log file

You might not realize it, but this is something that you now know how to do. Let's get started...

Rather than doing all of this in the terminal we are going to create a script file with all relevant commands. Use `vim` to create our new script file:

`$ vim generate_bad_reads_summary.sh`

We always want to start our scripts with a shebang line: 

`#!/bin/bash`

This line is the absolute path to the Bash interpreter. The shebang line ensures that the bash shell interprets the script even if it is executed using a different shell.

After the shebang line, we enter the commands we want to execute. First we want to move into our `~/ngs_course/dnaseq/data/untrimmed_fastq/` directory:

```
# enter directory with raw FASTQs
$ cd ~/ngs_course/dnaseq/data/untrimmed_fastq
```

And now we loop over all the FASTQs:

```
# loop over all FASTQ files
for filename in *.fastq.gz;
```

and we execute the commands for each loop:

```
do
  # tell us what file we're working on
  echo $filename;
  
  # grab all the bad read records into new file
  grep -B1 -A2 NNNNNNNNNN $filename > $filename-bad-reads.fastq;
``` 
  
We'll also count the number of these reads and put that in a new file, using the count (`-c`) flag of `grep`:

```
# grab the number of bad reads and write it to a summary file
grep -cH NNNNNNNNNN $filename >> bad-reads.count.summary;
done
```

If you've noticed, we slipped a new `grep` flag `-H` in there. This flag will report the filename along with the match string. This is useful for when we generate the summary file.



Exit out of `vim`, and voila! You now have a script you can use to assess the quality of all your new datasets. Your finished script, complete with comments, should look like the following:

```
#!/bin/bash 

# enter directory with raw FASTQs
cd ~/ngs_course/dnaseq/data/untrimmed_fastq

# loop over all FASTQ files
for filename in *.fastq.gz

do 
  echo $filename; 

  # grab all the bad read records
  grep -B1 -A2 NNNNNNNNNN $filename > $filename-bad-reads.fastq;

  # grab the number of bad reads and write it to a summary file
  grep -cH NNNNNNNNNN $filename >> bad-reads.count.summary;
done

```

To run this script, we simply enter the following command:

```
$ sh generate_bad_reads_summary.sh
```

To keep your data organized, let's move all of the bad read files and script out of our `untrimmed_fastq` directory into the `other` directory

`$ mkdir ~/ngs_course/dnaseq/other #in case you have not created this directory yet`

`$ mv ~/ngs_course/dnaseq/data/untrimmed_fastq/*bad* ~/ngs_course/dnaseq/other`



---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson was derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*
