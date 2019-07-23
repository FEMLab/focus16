from .run_all import get_and_check_ave_read_len_from_fastq
import os
import shutil
import logging as logger

def test_get_ave():
    reads = os.path.join(os.path.dirname(__file__), "test_data", "test_reads1.fq")
    test_result = get_and_check_ave_read_len_from_fastq(fastq1=reads, N=50, logger=logger)
    print(test_result)
    assert 150.0 == test_result
    return()
