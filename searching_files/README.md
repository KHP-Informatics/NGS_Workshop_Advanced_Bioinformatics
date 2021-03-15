---
title: "The Shell: Searching and Redirection"
author: "Alfredo Iacoangeli"
date: "15/03/2021"
---

## IMPORTANT: the paths to files and directories in this workshop are only examples. You will need use your own so please try not to copy and paste the commands but write them yourself matching the correct locations on your computer

Approximate time: 60 minutes

## Learning Objectives

* learn how to search for characters or patterns in a text file using the `grep` command
* learn about output redirection
* learn how to use the pipe (`|`) character to chain together commands

## Searching files

We went over how to search within a file using `less`. We can also
search within files without even opening them, using `grep`. Grep is a command-line
utility for searching plain-text data sets for lines matching a pattern or regular expression (regex).
Let's give it a try!

Suppose we want to see how many reads in our file `WES01_chr22m_R1.fastq.gz` are bad, with 6 consecutive Ns  
Let's search for the string NNNNNN: 

`$ cd ~/ngs_course/dnaseq/data/untrimmed_fastq`

`$ zcat WES01_chr22m_R1.fastq.gz > WES01_chr22m_R1.fastq #fastq.gz is a compressed format. I need to uncompress it first`

`$ grep NNNNNN WES01_chr22m_R1.fastq`

We get back a lot of lines.  What if we want to see the whole fastq record for each of these reads?
We can use the '-B' and '-A' arguments for grep to return the matched line plus one before (-B 1) and two
lines after (-A 2). Since each record is four lines and the second line is the sequence, this should
return the whole record.

`$ grep -B1 -A2 NNNNNN WES01_chr22m_R1.fastq`

for example:

	@HWI-ST330:304:H045HADXX:1:1101:1111:61397
	CACTTGTAAGGGCAGGCCCCCTTCACCCTCCCGCTCCTGGGGGANNNNNNNNNNANNNCGAGGCCCTGGGGTAGAGGGNNNNNNNNNNNNNNGATCTTGG
	+
	@?@DDDDDDHHH?GH:?FCBGGB@C?DBEGIIIIAEF;FCGGI#########################################################
	--
	@HWI-ST330:304:H045HADXX:1:1101:1106:89824
	CACAAATCGGCTCAGGAGGCTTGTAGAAAAGCTCAGCTTGACANNNNNNNNNNNNNNNNNGNGNACGAAACNNNNGNNNNNNNNNNNNNNNNNNNGTTGG
	+
	?@@DDDDDB1@?:E?;3A:1?9?E9?<?DGCDGBBDBF@;8DF#########################################################

****
**Exercise**

1. Search for the sequence TGCTTTCTTGA in `WES01_chr22m_R1.fastq`.
In addition to finding the sequence, have your search also return
the name of the sequence.
2. Search for that sequence in WES01_chr22m_R2.fastq.gz also.

****

## Redirection

We're excited we have all these sequences that we care about that we
just got from the FASTQ files. That is a really important motif
that is going to help us answer our important question. But all those
sequences just went whizzing by with grep. How can we capture them?

We can do that with something called "redirection". The idea is that
we're redirecting the output from the terminal (all the stuff that went
whizzing by) to something else. In this case, we want to print it
to a file, so that we can look at it later.

The redirection command for putting something in a file is `>`

Let's try it out and put all the sequences that contain 'NNNNNN'
from all the files in to another file called `bad_reads.txt`.

`$ grep -B1 -A2 NNNNNN WES01_chr22m_R1.fastq > bad_reads.txt`

`$ ls -l`

The prompt should sit there a little bit, and then it should look like nothing
happened. But type `ls`. You should have a new file called `bad_reads.txt`. Take
a look at it and see if it has what you think it should.

If we use '>>', it will append to rather than overwrite a file.  This can be useful for saving more than one search, for example:
 
`$ zcat WES01_chr22m_R2.fastq.gz > WES01_chr22m_R2.fastq `
 
`$ grep -B1 -A2 NNNNNN WES01_chr22m_R2.fastq >> bad_reads.txt`

`$ ls -l`

Since our `bad_reads.txt` file isn't a untrimmed_fastq file, we should move it to a different location within our directory. We decide to create a new folder called `other`, and move the `bad_reads.txt` to this `other` folder using the command `mv`. 

`$ mkdir ../other/`

`$ mv bad_reads.txt ../other/`

There's one more useful redirection command that we're going to show, and that's
called the pipe command, and it is `|`. It's probably not a key on
used very much on your keyboard. What `|` does is take the output that
scrolling by on the terminal and then can run it through another command.
When it was all whizzing by before, we wished we could just slow it down and
look at it, like we can with `less`. Well it turns out that we can! We pipe
the `grep` command to `less`

`$ grep -B1 -A2 NNNNNN WES01_chr22m_R1.fastq | less`

Now we can use the arrows to scroll up and down and use `q` to get out.

We can also do something tricky and use the command `wc`. `wc` stands for
*word count*. It counts the number of lines or characters. So, we can use
it to count the number of lines we're getting back from our `grep` command.
And that will magically tell us how many sequences we're finding.

`$ grep NNNNNN WES01_chr22m_R1.fastq | wc`

This command tells us the number of lines, words and characters in the file. If we
just want the number of lines, we can use the `-l` flag for `lines`.

`$ grep NNNNNN WES01_chr22m_R1.fastq | wc -l`

Redirecting is not super intuitive, but it's powerful for stringing
together these different commands, so you can do whatever you need to do.

The philosophy behind these commands is that none of them
really do anything all that impressive. BUT when you start chaining
them together, you can do some really powerful things 
efficiently. If you want to be proficient at using the shell, you must
learn to become proficient with the pipe and redirection operators:
`|`, `>`, `>>`.

## Practice with searching and redirection

Finally, let's use the new tools in our kit and a few new ones to examine our region annotation file, **chr22.genes.hg19.bed**, which we will be using later to find the genomic coordinates of all known exons on chromosome 22.

`$ cd ~/ngs_course/dnaseq/data`

Let's explore our `chr22.genes.hg19.bed` file a bit. What information does it contain?

`$ less chr22.genes.hg19.bed`

chr22   16279194        16279301        POTEH-chr22-16279195-16279301
chr22   16287253        16287937        POTEH-chr22-16287254-16287937
chr22   16448825        16449804        OR11H1-chr22-16448826-16449804
chr22   17071647        17073700        CCT8L2-chr22-17071648-17073700
chr22   17082800        17083105        psiTPTE22-chr22-17082801-17083105
chr22   17103730        17103787        psiTPTE22-chr22-17103731-17103787

> The BED file is a tab-delimited region annotation file often used in NGS analyses. For more information on this file format, check out the [Ensembl site](http://www.ensembl.org/info/website/upload/bed.html). 

The columns in the **BED file contain the genomic coordinates of a genomic region  (in our case exons) and then an optional number of optional tab-delimited fields. In our case the optional field contains the gene name and the exon genomic coordinates in a dash delimited format **. Note that sometimes an exon can be associated with multiple different genes and/or transcripts. 

Now that we know what type of information is inside of our gtf file, let's explore our commands to answer a simple question about our data: **how many total exons are present on chromosome 22 using `chr22.genes.hg19.bed`?**

To determine the number of total exons on chromosome 22, we are going to perform a series of steps:
	
	1. Subset the dataset to only include the genomic location information
	2. Extract only the genomic coordinates of RABL2B exons
	3. Remove duplicate regions
	4. Count the total number of exons
	
#### Subsetting dataset
We will define an exon by its genomic coordinates. Therefore, we only need the feature type and the genomic location (chr, start, stop, and strand) information to find the total number of exons. The columns corresponding to this information are 1, 3, 4, 5, and 7. 

'cut' is a program that will extract columns from files.  It is a very good command to know.  Let's first try out the 'cut' command on a small dataset (just the first 5 lines of chr22.genes.hg19.bed) to make sure we have the command correct:

`$ cut -f1,2,3 chr22.genes.hg19.bed | head -n 5`
   
'-f1,2,3' means to cut these fields (columns) from the dataset.  

chr22	16256331	16256677
chr22	16258184	16258303
chr22	16266928	16267095
chr22	16269872	16269943
chr22	16275206	16275277

The `cut` command assumes our data columns are separated by tabs (i.e. tab-delimited). The `chr22.genes.hg19.bed` is a tab-delimited file, so the default `cut` command works for us. However, data can be separated by other types of delimiters. Another common delimiter is the comma, which separates data in comma-separated value (csv) files. If your data is not tab delimited, there is a `cut` command argument (-d) to specify the delimiter.

Our output looks good, so let's cut these columns from the whole dataset (not just the first 5 lines) and save it as a file, **`chr22.genes.hg19.bed_cut`**:

`$ cut -f1,2,3 chr22.genes.hg19.bed > chr22.genes.hg19.bed_cut`

Check the cut file to make sure that it looks good using `less`. 

#### Extracting genomic coordinates of exon features
We only want the RABL2B exons, so let's use `grep` to only keep the exon lines and save to file, **`RABL2B_exons`**:

`$ grep RABL2B chr22.genes.hg19.bed > RABL2B_exons`

#### Removing duplicate exons
Now, we need to remove those exons that show up multiple times. As there are none in your bed file, let's create some duplicate lines:

`$ cat chr22.genes.hg19.bed | wc -l #this counts the number of lines in the file`

`$ cp chr22.genes.hg19.bed chr22.genes.hg19_plus5.bed #we make a copy of the file` 

`$ head -5 chr22.genes.hg19.bed >> chr22.genes.hg19_plus5.bed #this adds (appends) the first 5 lines of the file to itself`

`$ cat chr22.genes.hg19_plus5.bed | wc -l #the line count should be +5 now `

We can use a new tool, `sort`, to remove exons that show up more than once.  We can use the `sort` command with the `-u` option to return only unique lines.

`$ sort -u chr22.genes.hg19_plus5.bed | wc -l #the line count should be the original one now`

    



## Where can I learn more about the shell?

- Software Carpentry tutorial - [The Unix shell](http://software-carpentry.org/v4/shell/index.html)
- The shell handout - [Command Reference](http://files.fosswire.com/2007/08/fwunixref.pdf)
- [explainshell.com](http://explainshell.com)
- http://tldp.org/HOWTO/Bash-Prog-Intro-HOWTO.html
- man bash
- Google - if you don't know how to do something, try Googling it. Other people
have probably had the same question.
- Learn by doing. There's no real other way to learn this than by trying it
out.  Write your next paper in nano (really emacs or vi), open pdfs from
the command line, automate something you don't really need to automate.


**Commands, options, and keystrokes covered in this lesson**

```bash
grep
> (output redirection)
>> (output redirection, append)
| (output redirection, pipe)
wc
cut
sort
uniq (sort -u)
```

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*

* *The materials used in this lesson were derived from work that is Copyright © Data Carpentry (http://datacarpentry.org/). 
All Data Carpentry instructional material is made available under the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0).*

* *Adapted from the lesson by Tracy Teal. Contributors: Paul Wilson, Milad Fatenejad, Sasha Wood, and Radhika Khetani for Software Carpentry (http://software-carpentry.org/)*
