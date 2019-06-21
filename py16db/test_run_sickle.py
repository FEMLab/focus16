from run_sickle import run_sickle
import os
import shutil
from Bio import SeqIO

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def test_sickle_singles():
    sickle_test_dir = "tmp_sickle"
    fastq1 = "../riboSeed/tests/references/toy_reads1.fq"
    fastq2 = None
    if os.path.exists(sickle_test_dir):
        shutil.rmtree(sickle_test_dir)
    new_fastq1, new_fastq2 = run_sickle(fastq1=fastq1, fastq2 = None,
                                        output_dir=sickle_test_dir)
    assert new_fastq2 == None
    assert file_len(fastq1) > file_len(new_fastq1)
