from run_all import coverage
import os

def test_coverage():
    reads = "./plentyofbugs/test_data/test_reads1.fq"
    destination = "../test/0/"
    test_result = coverage(approx_length=2000000, fastq1=reads, SRA="SRR8443698", destination=destination, N=50)
    print(test_result)
    assert 0.735  == test_result
    return()
    
