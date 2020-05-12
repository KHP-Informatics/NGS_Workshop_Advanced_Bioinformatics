---
title: "Alignment with bwa and duplicate marking with picard"
author: "Alfredo Iacoangeli"
date: "16/03/2020"
---

## IMPORTANT: the paths to files and directories in this workshop are only examples. You will need use your own so please try not to copy and paste the commands but write them yourself matching the correct locations on your computer

## Learning Objectives:

* Understanding the alignment method BWA utilizes to align sequence reads to the reference genome
* Identifying the intricacies of alignment tools used in NGS analysis (parameters, usage, etc)
* Choosing appropriate BWA alignment parameters for our dataset
* Running BWA on multiple samples
* Marking duplicates with picard
* Filter BAM based on mapping quality and bitwise flags
* Some challenging exercises

## Read Alignment

Now that we have our quality-trimmed reads, we can move on to read alignment. We perform read alignment USING BWA MEM.

Please run bwa and spend a moment to see its options. Try to understand what they mean and research them on line if something is not clear.

```
$ bwa

Program: bwa (alignment via Burrows-Wheeler transformation)
Version: 0.7.17-r1188
Contact: Heng Li <lh3@sanger.ac.uk>

Usage:   bwa <command> [options]

Command: index         index sequences in the FASTA format
         mem           BWA-MEM algorithm
         fastmap       identify super-maximal exact matches
         pemerge       merge overlapping paired ends (EXPERIMENTAL)
         aln           gapped/ungapped alignment
         samse         generate alignment (single ended)
         sampe         generate alignment (paired ended)
         bwasw         BWA-SW for long queries

         shm           manage indices in shared memory
         fa2pac        convert FASTA to PAC format
         pac2bwt       generate BWT from PAC
         pac2bwtgen    alternative algorithm for generating BWT
         bwtupdate     update .bwt to the new format
         bwt2sa        generate SA from BWT and Occ

Note: To use BWA, you need to first index the genome with `bwa index'.
      There are three alignment algorithms in BWA: `mem', `bwasw', and
      `aln/samse/sampe'. If you are not sure which to use, try `bwa mem'
      first. Please `man ./bwa.1' for the manual.


```

## BWA MEM 

BWA (please visit the developer's github https://github.com/lh3/bwa) is a software package for mapping DNA sequences against a large reference genome, such as the human genome. It consists of three algorithms: BWA-backtrack, BWA-SW and BWA-MEM. The first algorithm is designed for Illumina sequence reads up to 100bp, while the rest two for longer sequences ranged from 70bp to a few megabases. BWA-MEM and BWA-SW share similar features such as the support of long reads and chimeric alignment, but BWA-MEM, which is the latest, is generally recommended as it is faster and more accurate. BWA-MEM also has better performance than BWA-backtrack for 70-100bp Illumina reads.

For all the algorithms, BWA first needs to construct the FM-index for the reference genome (the index command). Alignment algorithms are invoked with different sub-commands: aln/samse/sampe for BWA-backtrack, bwasw for BWA-SW and mem for the BWA-MEM algorithm.

Please run bwa and spend a moment to see its options. Try to understand what they mean and research them on line if something is not clear.

```

$ bwa mem

Algorithm options:

       -t INT        number of threads [1]
       -k INT        minimum seed length [19]
       -w INT        band width for banded alignment [100]
       -d INT        off-diagonal X-dropoff [100]
       -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [1.5]
       -y INT        seed occurrence for the 3rd round seeding [20]
       -c INT        skip seeds with more than INT occurrences [500]
       -D FLOAT      drop chains shorter than FLOAT fraction of the longest overlapping chain [0.50]
       -W INT        discard a chain if seeded bases shorter than INT [0]
       -m INT        perform at most INT rounds of mate rescues for each read [50]
       -S            skip mate rescue
       -P            skip pairing; mate rescue performed unless -S also in use

Scoring options:

       -A INT        score for a sequence match, which scales options -TdBOELU unless overridden [1]
       -B INT        penalty for a mismatch [4]
       -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
       -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
       -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
       -U INT        penalty for an unpaired read pair [17]

       -x STR        read type. Setting -x changes multiple parameters unless overridden [null]
                     pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0  (PacBio reads to ref)
                     ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0  (Oxford Nanopore 2D-reads to ref)
                     intractg: -B9 -O16 -L5  (intra-species contigs to ref)

Input/output options:

       -p            smart pairing (ignoring in2.fq)
       -R STR        read group header line such as '@RG\tID:foo\tSM:bar' [null]
       -H STR/FILE   insert STR to header if it starts with @; or insert lines in FILE [null]
       -o FILE       sam file to output results to [stdout]
       -j            treat ALT contigs as part of the primary assembly (i.e. ignore <idxbase>.alt file)
       -5            for split alignment, take the alignment with the smallest coordinate as primary
       -q            don't modify mapQ of supplementary alignments
       -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility) []

       -v INT        verbosity level: 1=error, 2=warning, 3=message, 4+=debugging [3]
       -T INT        minimum score to output [30]
       -h INT[,INT]  if there are <INT hits with score >80% of the max score, output all in XA [5,200]
       -a            output all alignments for SE or unpaired PE
       -C            append FASTA/FASTQ comment to SAM output
       -V            output the reference FASTA header in the XR tag
       -Y            use soft clipping for supplementary alignments
       -M            mark shorter split hits as secondary

       -I FLOAT[,FLOAT[,INT[,INT]]]
                     specify the mean, standard deviation (10% of the mean if absent), max
                     (4 sigma from the mean if absent) and min of the insert size distribution.
                     FR orientation only. [inferred]

Note: Please read the man page for detailed description of the command line and options.
```

### Making the index

When building the reference index, BWA will generate several files. Let's create a folder for the reference and its index files and then run bwa to generate the index files.

```
$ mkdir -p ~/ngs_course/dnaseq/data/reference

$ mv ~/ngs_course/dnaseq/data/hg19.fa.gz ~/ngs_course/dnaseq/data/reference/

$ bwa index ~/ngs_course/dnaseq/data/reference/hg19.fa.gz
```

Indexing the reference will take a while (~45 minutes) so I would suggest that after running it, you open another terminal and ssh onto the OpenStack instance so that you can keep working on the new terminal, of course you have to wait for the index to be ready if you want to go ahead with the bwa alignment. BWA index will generate a number of files in ~/ngs_course/dnaseq/data/reference/. Please list them and google what they are.

```
$ ls ~/ngs_course/dnaseq/data/reference

hg19.fa.gz  hg19.fa.gz.amb  hg19.fa.gz.ann  hg19.fa.gz.bwt  hg19.fa.gz.pac  hg19.fa.gz.sa

```

### Read group infos

We will use the following read group infos for the alignment, they are the same we used for the Galaxy workshop:


-Read group identifier (ID): HWI-D0011.50.H7AP8ADXX.1.WES01
-Read group identifier (SM): WES01
-Read group identifier (PL): ILLUMINA
-Read group identifier (LB): nextera-wes01-blood
-Read group identifier (PU): HWI-D00119
-Read group identifier (DT): 2017-02-23

If you want to read more about RG infos please visit https://galaxyproject.org/learn/galaxy-ngs101/#read-groups


### Running BWA mem with the RG information

```
$ mkdir ~/ngs_course/dnaseq/data/aligned_data

$ bwa mem -t 4 -v 1 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50  ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq/data/trimmed_fastq/WES01_chr22m_trimmed_R_1P.fastq ~/ngs_course/dnaseq/data/trimmed_fastq/WES01_chr22m_trimmed_R_2P.fastq > ~/ngs_course/dnaseq/data/aligned_data/WES01_chr22m.sam

```

You should have a directory tree setup similar to that shown below. it is best practice to have all files you intend on using for your workflow present within the same directory. In our case, we have our original FASTQ files and post-trimming data generated in the previous section. We also have all reference data files that will be used in downstream analyses.

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
	│       └── WES01_chr22m.sam
	├── meta
	├── results
	└── logs
```

Change directories into the `aligned_data` folder. 

```
$ cd ~/ngs_course/dnaseq/data/aligned_data
```

Now we need to convert the sam file into bam format, sort it and generate an index using samtools.
```
$ samtools view -h -b WES01_chr22m.sam > WES01_chr22m.bam #IMPORTANT: depending on your samtools version, you might need to add the -S option to make samtools accept the sam input file 

$ samtools sort WES01_chr22m.bam > WES01_chr22m_sorted.bam

$ samtools index WES01_chr22m_sorted.bam #This will generate a .bai index file

$ ls

WES01_chr22m.bam  WES01_chr22m.sam  WES01_chr22m.bam  WES01_chr22m_sorted.bam.bai
```

## Post Alignment QC and Filtering

We will Use 1) NGS: Picard > MarkDuplicates and then 2) SAMtools > Filter SAM or BAM, output SAM or BAM files to make a new BAM file which flags duplicate molecules and excludes poor quality reads.

### MarkDuplicates

We will use Picard to mark duplicated reads. This tool examines aligned records in the supplied SAM or BAM dataset to locate duplicate molecules. All records are then written to the output file with the duplicate records flagged. Two files are produced, 1) the new BAM file with duplicate reads marked and 2) a metrics file summarising the number of duplicate reads found.

If you want to read more about pcr duplicates visit http://www.cureffi.org/2012/12/11/how-pcr-duplicates-arise-in-next-generation-sequencing/

```
$ picard MarkDuplicates I=WES01_chr22m_sorted.bam O=WES01_chr22m_sorted_marked.bam M=marked_dup_metrics.txt

$ samtools index WES01_chr22m_sorted_marked.bam
```


### Filter BAM based on mapping quality and bitwise flags using samtools

We are going to filter the reads according to the following criteria:

-Set Minimum MAPQ quality score : 20
-Filter on bitwise flag: yes a. Skip alignments with any of these flag bits set
  i. The read is unmapped
  ii. The alignment or this read is not primary
  iii. The read fails platform/vendor quality checks
  iv. The read is a PCR or optical duplicate

To interpret the samtools flags we are using to filter please see https://broadinstitute.github.io/picard/explain-flags.html


```
$ samtools view -F 1796  -q 20 -o WES01_chr22m_sorted_filtered.bam WES01_chr22m_sorted_marked.bam

$ samtools index WES01_chr22m_sorted_filtered.bam

```


## Excercise

1)So far we have followed the main steps you performed with Galaxy during the first module. We have not gone through a few of the alignment statistics part yet and we will not provide detailed instructions about how to do it in the command line.
As an exercise please try to perform the Alignment Statistics steps of the Galaxy workshops on your terminal by yourself. You alignment statistics analysis should include the following steps:
 
 - Flagstats
 - Viewing the BAM File
 - Generate alignment statistics per chromosome
 - Determine the distribution of insert sizes
 - Depth of Coverage

The first three steps can be performed with samtools. The fourth( Determine the distribution of insert sizes) will require you use picard. Picard is less intuitive to use comparing with samtools and I would suggest you look for its documentation on line. The last one (Depth of Coverage) will require you to use bedtools.
You are welcome to go though the Galaxy workshop material to refresh your memory. If you don't have the material anymore please contact Alfredo on Slack or Email (alfredo.iacoangeli@kcl.ac.uk).

2) Take a moment to update the README for the dnaseq folder (hint: use vim, or nano or any text editor of your choice to create the file). Give a short update of the project and brief descriptions of the types of file you have generated within each of the sub-directories. Please take note of the current size of the project. The total storage available on your virtual machine if 40Gigabytes

---
*This lesson has been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*






