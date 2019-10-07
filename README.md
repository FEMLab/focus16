[![Build Status](https://travis-ci.com/FEMLab/focusdb.svg?branch=master)](https://travis-ci.com/FEMLab/focusdb)[![PyPI version](https://badge.fury.io/py/focusDB.svg)](https://badge.fury.io/py/focusDB)[![DOI](https://zenodo.org/badge/191556439.svg)](https://zenodo.org/badge/latestdoi/191556439)


# focusDB
## High resolution 16S database construction from correctly assembled rDNA operons

## Description
focusDB is a package built for the construction of species-specific, high-resolution 16S rDNA databases.
It does so with through the use of riboSeed, a pipeline for the use of ribosomal flanking regions to improve bacterial genome assembly.

## What it does:
focusDB sets up a directory in your home folder to store  downloaded SRAs. This makes rerunning easier, as all the raw data is centralized for reuse. This keeps the data transfer load to a minimum.  A manifest file keeps track of the runs already downloaded.

Given an organism, focusDB does the following:
	- Identify all the whole-genome seqeuncing SRAs available for that species
	- Download all potential complete reference genomes.
	For each SRA:
		- use `plentyofbugs` to identify which reference genome would be the closest
		- run QC on reads, downsampling if neccessary
		- run riboSeed assembly
	Extract all the 16S seqeunces from reassembly


## Installation
#### Installing focusDB
focusDB available via pypi. We recommend installing within a python environment.

```
conda create --name focusDBenv python=3.5 seqtk sickle-trim sra-tools riboseed mash skesa barrnap parallel-fastq-dump iqtree mafft
conda activate focusDBenv
pip install focusDB

```

#### Packages needed for additional features:
Optionally, to use the trimming alignment feature, TrimAl must be installed from github https://github.com/scapella/trimal.  For re-generating test data, ART read simulator must also be installed.


## Usage
####  reassemble SRAs and extract potentially novel 16S sequences
This will go through the process of setting up the `.focusDB` dir (defaults to your home directory), downloading up to 30 reference complete genomes, downloading up to 5 WGS SRAs, finding the closes referece for each of the 5 SRAs, assembling, and extracting the 16S sequences.  Say we have a 4 core computer with 16 gb ram, we spilt it so the assemblies each use half the resources.
```
focusDB --output_dir ./focusdb_ecoli/ --n_SRAs 5 --n_references 30 --memory 8 --cores 2 --njobs 2 --organism_name "Escherichia coli"
```
#### Optional downstream analyses
```
# build E. coli specific DB from E colis in Silva and our new seqeunces
combine-focusdb-and-silva  -d ~/Downloads/SILVA_132_SSUParc_tax_silva.fasta  -o ecolidb.fasta  -n "Escherichia coli" -S ./focusdb_ecoli/ribo16s.fasta
# Align sequences and trim  the alignment
align-and-trim-focusdb -i ecolidb.fasta --out_prefix aligned_ecolidb
# Calculate the per-column shannon entropy of the trimmed alignment.
calculate-shannon-entropy calculate-shannon-entropy.py -i aligned_ecolidb.mafft.trimmed > ecoli_entropy
```



##### `focusDB`
###### Required Arguments
```
[--organism_name]: The species of interest, input within quotes.
[--n_references]: Maximum number of reference genomes the user wishes to download.
[--n_SRAs]: Maximum number of SRAs the user wishes to download.
```
###### Optional Arguments:
```
[--sra_list]: Uses a user-given list of SRA accessions instead of obtaining SRA accessions from the pipeline.
[--version]: Returns focusDB version number.
[--approx_length]: Uses a user-given genome length as opposed to using reference genome length.
[--sraFind_path]: Path to pre-downloaded sraFind-All-biosample-with-SRA-hits.txt file.
[--prokaryotes]: Path to pre-downloaded prokaryotes.txt file.
[--get_all]: If one SRA has two runs, downloads both.
[--cores]: The number of cores the user would like to use for focusDB. Specifically, riboSeed and plentyofbugs can be optimized for thread usage.
[--memory]: As with [--cores], RAM can be optimized for focusDB.
[--maxcov]: The maximum read coverage for SRA assembly. Downsamples to this coverage if the coverage exceeds it.
[--example_reads]: Input of user-given reads.
[--subassembler]: Choice of mash or skesa for subassembly in riboSeed.
```

### Included Utilities:
#### `combine-focusdb-and-silva`
Use this script to combine silva  and focusDB seqeunces for a given organism name.
#### `align-and-trim-focusdb`
```
usage: align-and-trim-focusdb [-h] -i INPUT -o OUT_PREFIX

Given a multiple sequence database (from combine-focusdb-and-silva, generate
an alignment with mafft and trim to median sequence. Requires mafft and TrimAl

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        multifasta input file
  -o OUT_PREFIX, --out_prefix OUT_PREFIX
                        prefix for msa and trimmed msa
```

#### `calculate-shannon-entropy`
```
usage: calculate-shannon-entropy [-h] -i INPUT

Given a trimmed multiple sequence alignment (from align-and-trim-focusdb,
calculate shannon entropy

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        trimmed MSA
```


## Beast Mode
Reassembling hundreds or thousands of genomes can eat into a ton of resources, so we try to make this a feasible. For optimal running
1) install Aspera connect.  NCBI uses this to transfer data faster than the default connections with http or ftp.
2) set up your cache somewhere with lots of space and fast I/O.  run `vdb-config -i`, following the instructions here: https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration
3) Use the focusDB-prefetch command to get your data. This will download the data as `.sra`, and add them to your cache.  This can be configured to run the requests in batches.
4) Run focusDB: the calls to `fasterq-dump` will now go to the local copy of the `.sra` in the cache area.


## Test Data
### Unit tests
Testing is done with the `nose` package. Generate the test data with
```
nosetests  pyfocusDB/generator.py
```
and run the unit tests with

```
nosetests pyfocusDB/ -v
```

Note  that `generator.py` requires ART to generate synthetic.
{https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm}

### Running on test datasets



## Bugs

### Numpy
If you get a failure running riboSeed about `dependencies not installed:["numpy"]`, try running `python -c "import numpy as np"`. If you get an error about multiple versions being installed, repeat `pip uninstall numpy` repeatedly until no more versions remain in your environment.  Then, reinstall numpy and try again.

### OpenBlas on MacOS
If you get an error about openblas, try upgrading the one chosen by conda with:
```
conda install openblas=0.2.19
```


### Fuzzy matching organisms for `plentyofbugs`
The default behavior for identifying organisms of interest from NCBI's `prokaryotes.txt` is to find lines starting with `--organism_name`.  This is intentional, as the names are poorly defined in the file, and  this allows us to capture a whole genus if desired.  We have never come across a case where this happens, but this could have the consequence of including undesired organisms, if they start with the same characters. For instance, an `--organism_name` of `dog` would also match `dogfish`.  If you notice undesired organisms included in the log file in our output directory,  you will have to manually select the lines of interest from `prokaryotes.txt`, save that as an alterative file, and use the `--prokaryotes` argument to provide this edited version to focusDB.
