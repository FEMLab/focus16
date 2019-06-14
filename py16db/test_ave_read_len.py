from run_all import get_ave_read_len_from_fastq
import os
import shutil

def test_get_ave():
    reads = "./plentyofbugs/test_data/test_reads1.fq"
    test_result = get_ave_read_len_from_fastq(fastq1=reads, N=50)
    assert ["read length is OK : 150.0 bp"] == test_result
    return()
