from run_all import get_coverage, downsample
import os
import shutil
from nose.tools.nontrivial import with_setup


def test_coverage():
#genome is built from NC_0176759.1 ~150000bp
    reads = os.path.join(os.path.dirname(__file__), "test_data", "test_reads1.fq")
    test_result = get_coverage(approx_length=150000, fastq1=reads)
    print(test_result)
    assert 9.8  == test_result
    return()

def setup_func():
    "set up test fixtures"
    if os.path.exists("./py16db/downsample_test_result/"):
        shutil.rmtree("./py16db/downsample_test_result/")
    pass


def teardown_func():
    "tear down test fixtures"
    shutil.rmtree("./py16db/downsample_test_result/")

@with_setup(setup_func, teardown_func)
def test_downsample():
#genome is built from NC_0176759.1 ~150000bp
    reads = os.path.join(os.path.dirname(__file__), "test_data", "test_reads1.fq")
    os.makedirs("./py16db/downsample_test_result/")
    os.makedirs("./py16db/downsample_test_result/downsampled/")
    destination = "./py16db/downsample_test_result/"
    reads1, reads2 = downsample(approx_length=150000, fastq1=reads, fastq2=reads,
                                destination=destination, maxcoverage=2)
    down_cov = get_coverage(approx_length=150000, fastq1=reads1)
    assert 1.993 == down_cov
    return()
    
