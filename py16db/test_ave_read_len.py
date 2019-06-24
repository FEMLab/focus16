from run_all import get_ave_read_len_from_fastq
import os
import shutil


def test_get_ave():
    reads = os.path.join(os.path.dirname(__file__),"data", "test_data", "test_reads1.fq")
    test_result = get_ave_read_len_from_fastq(fastq1=reads, N=50)
    print(test_result)
    assert 150.0 == test_result
    return()
