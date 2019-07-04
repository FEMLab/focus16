# 16db 
## High resolution 16S database construction from correctly assembled rDNA operons

## Description
16db is a package built for the construction of species-specific, high-resolution 16S rDNA databases. 
It does so with through the use of riboSeed, a pipeline for the use of ribosomal flanking regions to improve bacterial genome assembly.
riboSeed allows the correct assembly of multiple rDNA operons within a single genome. 16db uses various tools around and including
riboSeed to take an input of arguments listed below and produces a file containing all 16S sequences from draft full genomes available for that species.  


## Installation
Packages required for 16db:  
```conda install seqtk sickle-trim sra-tools riboseed mash skesa plentyofbugs barrnap```
```pip install riboseed```
```pip install plentyofbugs```

## Arguments

###### Required: 
'[-o]' output directory  
'[-n]' organism name in quotes
'[-l]' approx length of organism genome
'[-g]' directory containing genomes in fasta/ or desired directory
'[-S]' number of strains for comparison to find reference genome

###### Optional:
'[-s]' path to sraFind-CompleteGenome-biosample-with-SRA-hits.txt 
'[--single_SRA]' run pipeline for one given SRA only
'[-p]' path to prokaryotes.txt 
'[--get_all]' get both SRAs if organism has two
'[--cores]' the amount of cores your computer has/ you want to use
'[--maxcov]' maximum coverage of reads, downsamples if over
'[--example_reads]' use own reads: forwardreads.fq reversereads.fq OR reads.fq 