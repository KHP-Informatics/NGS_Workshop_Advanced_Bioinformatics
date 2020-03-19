---
title: "Performing variant calling with Freebayes and variant filtering with vcflib"
author: "Alfredo Iacoangeli"
date: "16/03/2020"
---

## IMPORTANT: the paths to files and directories in this workshop are only examples. You will need use your own so please try not to copy and paste the commands but write them yourself matching the correct locations on your computer


Accurate and consistent variant calling requires statistical modelling and is essential for the clinical implementation of NGS. However, many programs are available for calling variants and their concordance varies. Furthermore, variants have different levels of confidence due to differences in data quality. For variants with intermediate confidence levels, it is difficult to separate true variation from artefacts that arise from many factors such as sequencing error, misalignment and inaccurate base quality scores. As a result, the evidence for variant calls requires scrutiny and caution should be used when interpreting positive and negative findings especially for indels which are more error prone. At the end of this exercise you will be able to:

-Use a a state of the art software (Freebayes) to call small variants (SNVs and indels)
-Describe the contents of variant call files
-Generate and interpret basic variant quality control parameters
-Use quality control filters to exclude or flag variants with low confidence
-See: https://wiki.galaxyproject.org/Learn/GalaxyNGS101#Finding_variants

## Variant calling with Freebayes


FreeBayes is a Bayesian genetic variant detector designed to find small polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels (insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex events (composite insertion and substitution events) smaller than the length of a short-read sequencing alignment. See https://github.com/ekg/freebayes for details on FreeBayes.
