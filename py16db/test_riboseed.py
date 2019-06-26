from run_all import run_riboseed
import os
import shutil
from nose.tools.nontrivial import with_setup


def setup_func():
    "set up test fixtures"
    if os.path.exists("./py16db/riboSeed"):
        shutil.rmtree("./py16db/riboSeed/")
    pass

def teardown_func():
    "tear down test fixtures"
    shutil.rmtree("./py16db/riboSeed/")


@with_setup(setup_func, teardown_func)

def test_riboseed():
    readsf="./py16db/example_data/reads1.fq"
    readsr="./py16db/example_data/reads2.fq"
    output_dir="py16db/riboSeed/"
    os.makedir=(output_dir)
    sra="./py16db/example_data/NC_013928.1.fna"
          
    run_riboseed(sra=sra, readsf=readsf,
                 readsr=readsr, cores="4", threads="1",
                 output=output_dir)
    
    
    assert os.path.exists("./py16db/riboSeed/seed/final_long_reads/")

    return()
