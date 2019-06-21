from run_all import get_coverage, downsample
import os


def test_coverage():
#genome is built from NC_0176759.1 ~150000bp
    reads = "./plentyofbugs/test_data/test_reads1.fq"
    test_result = get_coverage(approx_length=150000, fastq1=reads)
    print(test_result)
    assert 9.8  == test_result
    return()

def test_downsample():
#genome is built from NC_0176759.1 ~150000bp
    reads = "./plentyofbugs/test_data/test_reads1.fq"
    destination = "../test/0/"
    reads1, reads2 = downsample(approx_length=150000, fastq1=reads, fastq2=reads,
                             destination=destination, maxcoverage=2)
    down_cov = get_coverage(approx_length=150000, fastq1=reads1)
    assert 1.993 == down_cov
    return()
    
