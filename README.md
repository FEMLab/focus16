# 16db
## High resolution 16S database construction from correctly assembled rDNA operons

## Description
16db is a package built for the construction of species-specific, high-resolution 16S rDNA databases.
It does so with through the use of riboSeed, a pipeline for the use of ribosomal flanking regions to improve bacterial genome assembly.
riboSeed allows the correct assembly of multiple rDNA operons within a single genome. 16db uses various tools around and including
riboSeed to take an input of arguments listed below and produces a file containing all 16S sequences from draft full genomes available for that species.


## Installation
###### Installing 16db
TODO: get pip install worrking
```
pip install 16db
```

###### Packages required for 16db:
```
conda install python=3.5 seqtk sickle-trim sra-tools riboseed mash skesa barrnap parallel-fastq-dump iqtree
```
Optionally, to use the trimming alignment feature, TrimAl must be installed from github https://github.com/scapella/trimal.  For re-generating test data, ART read simulator must also be installed.


## Usage
###### Example
```
# reassemble SRAs and extract potentially novel 16S sequnces
16db --output_dir ./focusdb_ecoli/ -g ./escherichia_genomes/ --n_SRAs 5 --n_references 30 --memory 8 --cores 4 --organism_name "Escherichia coli"
## Optional downstream analyses
# build E. coli specific DB from E colis in Silva and our new seqeunces
combine-focusdb-and-silva  -d ~/Downloads/SILVA_132_SSUParc_tax_silva.fasta  -o ecolidb.fasta  -n "Escherichia coli" -S ./focusdb_ecoli/ribo16s.fasta
# Align sequences and trim  the alignment
align-and-trim-focusdb -i ecolidb.fasta --out_prefix aligned_ecolidb
# Calculate the per-column shannon entropy of the trimmed allignment.
calculate-shannon-entropy calculate-shannon-entropy.py -i aligned_ecolidb.mafft.trimmed > ecoli_entropy
```



##### `16db`
This will go through the process of getting the list of assemblies that are associated with SRAs, downloading up to 5 SRAs,  finding the closes referece for each of the 5 SRAs, assembling, and extracting the 16S sequences.



###### Required Arguments
```
[--organism_name]: The species of interest, input within quotes.
[--nstrains]: The number of reference genomes and the number of SRAs the user wishes to download.
[--output_dir]: The output directory.
[--genomes_dir]: The output directory for which to store reference genomes, or a preexisting directory containing genomes the user wishes to use as reference genomes.
```
###### Optional Arguments:
```
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
```

### Included Utilities:
#### `combine-focusdb-and-silva`
Use this script to combine silva  and 16db seqeunces for a given organism name.
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

## Test Data
### Unit tests
Testing is done with the `nose` package. Generate the test data with
```
nosetests  py16db/generator.py
```
and run the unit tests with

```
nosetests py16db/ -v
```

Note  that `generator.py` requires ART to generate synthetic.
{https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm}

### Running on test datasets




## Bugs

### OpenBlas on MacOS
If you get a failure running riboSeed about `dependencies not installed:["numpy"]`, try running `python -c "import numpy as np"`. If you get an error about openblas, try upgrading the one chosen by conda with:
```
conda install openblas=0.2.19
```


### Fuzzy matching organisms for `plentyofbugs`
The default behavior for identifying organisms of interest from NCBI's `prokaryotes.txt` is to find lines starting with `--organism_name`.  This is intentional, as the names are poorly defined in the file, and  this allows us to capture a whole genus if desired.  We have never come across a case where this happens, but this could have the consequence of including undesired organisms, if they start with the same characters. For instance, an `--organism_name` of `dog` would also match `dogfish`.  If you notice undesired organisms included in the log file in our output directory,  you will have to manually select the lines of interest from `prokaryotes.txt`, save that as an alterative file, and use the `--prokaryotes` argument to provide this edited version to 16db.
