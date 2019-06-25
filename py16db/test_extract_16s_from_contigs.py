from run_all import extract_16s_from_contigs
from nose.tools.nontrivial import with_setup
import os
import shutil

#def setup_func():
#    "set up test fixtures"
#    if os.path.exists("./py16db/extract_test/"):
#       shutil.rmtree("./py16db/extract_test/")
#    pass
#
#def teardown_func():
#    "tear down test fixtures"
#    shutil.rmtree("./py16db/extract_test/")
#
#@with_setup(setup_func, teardown_func)
def test_extract():
    os.makedirs("./py16db/extract_test/")
    input_contigs = os.path.join(os.path.dirname(__file__),  "data", "test_data", "riboseed", "genomes", "NC_018089.1.fna")
    barr_out="./py16db/extract_test/barrnap"
    output="./py16db/extract_test/ribo16"
    test_result=extract_16s_from_contigs(input_contigs=input_contigs,
                                         barr_out=barr_out, output=output)
    assert []  == test_result



    

