# D2R_codes

This is the official repository for D2R, a tool that measures how repetitive a DNA sequence is. You give it a DNA sequence (or any sequence made of a limited set of characters), and it gives you a score — the higher the score, the more repetitive the sequence. You can also slide a window along a longer sequence, like an entire genome, to see how repetitiveness changes at different positions.

We currently offer several ways to use D2R, and we welcome contributions that make it even easier to access and run.

## 1. D2RReadFilter
A java package for filtering of sequencing reads. 
Only reads with D2R value higher than the threshold will be kept.

* Usage: 
```
# pt pc pg pa are the frequency of the four bases in your genome
java -jar D2RReadFilter    \path\to\read\fasta\file    kmersize    threshold    pt pc pg pa
```

Example:
```
java -jar D2RReadFilter.jar    reads.fa    7    0.00015    17 34 34 17
```

## 2. D2RGenomeScanner
A genome scanner spotting highly repetitive regions on genome.

The program is written in C++ with Visual Studio 2017.
Only standard libraries are used in the single cpp code file.

You can remove the #include "stdafx.h" line in the front of `D2R_codes/D2RGenomeScanner/D2Rscanner/D2Rscanner.cpp`
to compile it with other compilers or IDEs.

* Usage:
Build it as an executable file `D2Rscanner` and run it
```
D2Rscanner \path\to\genome\fasta\file  window_size  kmer_size
```

Example:
```
D2Rscanner E:\code\D2R\CRISPR\genomes\NC_010162.fa  1000  7
```


## 3. D2R python package
A python library for hacking. Written in swig, supporting Win64 and Linux64 environment currently.

Copy files under `release_win_x64` or `release_linux_x64` to your working directory and 
`import pyd2r`, D2R function is available with `def D2R_nNorm(seq, k, pt, pc, pg, pa)`
