from .run_all import alignment
from nose.tools.nontrivial import with_setup
import os
import shutil

def setup_func():
   "set up test fixtures"
   if os.path.exists("./py16db/alignment_test/"):
       shutil.rmtree("./py16db/alignment_test/")
       pass

def teardown_func():
    "tear down test fixtures"
    shutil.rmtree("./py16db/alignment_test/")

@with_setup(setup_func, teardown_func)
def test_alignment():
    test_output = "./py16db/alignment_test/"
    os.makedirs(test_output)
    test_fasta = "./py16db/test_data/test_16s_multiline.fasta"
    test_result = alignment(fasta=test_fasta, output=test_output)
    print(test_result)
    assert os.path.exists("./py16db/alignment_test/alignment/MSA.fasta")  
