---
title: "Automated bash pipeline"
author: "Alfredo Iacoangeli"
date: "16/03/2020"
---

## Learning Objectives:
* Learn how to automate your pipeline using a bash script 
* Use a `for loop` to run the pipeline on multiple files



Suppose you want to run a few different commands in succession – this is called “automating tasks”. See http://xkcd.com/1205/

Key caveat: none of the commands require any user live input (‘yes’, ‘no’, etc)

1) Figure out what commands you want to run.
2) Put them in a file using a text editor, like nano or pico, or TextEdit, or TextWrangler, or vi, or emacs.
3) Run the file.
4) For example, try putting the following shell commands in a bash script called ‘test.sh’:
