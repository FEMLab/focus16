from .run_all import pob
import os
import shutil
from nose.tools.nontrivial import with_setup


def setup_func():
    "set up test fixtures"
    if os.path.exists("./py16db/pob_test_result/"):
        shutil.rmtree("./py16db/pob_test_result/")
    pass

def teardown_func():
    "tear down test fixtures"
    shutil.rmtree("./py16db/pob_test_result/")


@with_setup(setup_func, teardown_func)
def test_pob():
    plasmids = os.path.join(os.path.dirname(__file__), "test_data", "plasmids", "")

    reads = os.path.join(os.path.dirname(__file__), "test_data", "test_reads1.fq")
    
    os.makedirs("./py16db/pob_test_result/")
    output_dir="./py16db/pob_test_result/plentyofbugs"
   
    test_result = pob(genomes_dir=plasmids, readsf=reads, output_dir=output_dir)
    print(test_result)
    assert ["./data/test_data/plasmids/NC_009837.1.fasta", "0.0340249"] == test_result
    return()


    
