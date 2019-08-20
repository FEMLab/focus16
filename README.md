# 16db 
## High resolution 16S database construction from correctly assembled rDNA operons

## Description
16db is a package built for the construction of species-specific, high-resolution 16S rDNA databases. 
It does so with through the use of riboSeed, a pipeline for the use of ribosomal flanking regions to improve bacterial genome assembly.
riboSeed allows the correct assembly of multiple rDNA operons within a single genome. 16db uses various tools around and including
riboSeed to take an input of arguments listed below and produces a file containing all 16S sequences from draft full genomes available for that species.  


## Installation
###### Installing 16db
```pip install 16db```

###### Packages required for 16db:  
```conda install seqtk sickle-trim sra-tools riboseed mash skesa barrnap prallel-fastq-dump```
```git clone https://github.com/nickp60/open_utils.git``` #For extractRegion command
```pip install plentyofbugs```



## Usage
###### Arguments
###### Required: 
[--organism_name]: The species of interest, input within quotes.  
[--nstrains]: The number of reference genomes and the number of SRAs the user wishes to download. 
[--output_dir]: The output directory. 
[--genomes_dir]: The output directory for which to store reference genomes, or a preexisting directory containing genomes the user wishes to use as reference genomes. 

###### Optional:
[--sra_list]: Uses a user-given list of SRA accessions instead of obtaining SRA accessions from the pipeline. 
[--version]: Returns 16db version number.
[--approx_length]: Uses a user-given genome length as opposed to using reference genome length. 
[--sraFind_path]: Path to pre-downloaded sraFind-All-biosample-with-SRA-hits.txt file. 
[--prokaryotes]: Path to pre-downloaded prokaryotes.txt file.
[--get_all]: If one SRA has two accessions, downloads both. 
[--cores]: The number of cores the user would like to use for 16db. Specifically, riboSeed and plentyofbugs can be optimized for thread usage.
[--memory]: As with [--cores], RAM can be optimized for 16db. 
[--maxcov]: The maximum read coverage for SRA assembly. Downsamples to this coverage if the coverage exceeds it. 
[--example_reads]: Input of user-given reads.
[--subassembler]: Choice of mash or skesa for subassembly in riboSeed.

###### Example
``` 16db -o ~/Desktop/escherichia/ -g ~/Desktop/escherichia/genomes/ -S 10 --memory 8 --cores 4 -n "Escherichia coli" ```

## Test Data
16db/py16db/generator.py is a script that generates the test data used for the
various tests in this repo

Required:
ART - generates synthetic next-gen reads.
{https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm}
Also move art_bin_MountRainier/ into 16db/py16db/




## Description
16db is a package built for the construction of species-specific, high-resolution 16S rDNA databases. 
It does so with through the use of riboSeed, a pipeline for the use of ribosomal flanking regions to improve bacterial genome assembly.
riboSeed allows the correct assembly of multiple rDNA operons within a single genome. 16db uses various tools around and including
riboSeed to take an input of arguments listed below and produces a file containing all 16S sequences from draft full genomes available for that species.  


## Installation
Packages required for 16db:  
```conda install seqtk sickle-trim sra-tools riboseed mash skesa barrnap prallel-fastq-dump```
```git clone https://github.com/nickp60/open_utils.git``` #For extractRegion command
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

## Test Data
16db/py16db/generator.py is a script that generates the test data used for the
various tests in this repo

Required:
ART - generates synthetic next-gen reads.
{https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm}
Also move art_bin_MountRainier/ into 16db/py16db/

Repo containing test data, we just need the plentyofbugs/test_data/plasmids folder.
{https://github.com/nickp60/plentyofbugs.git}

