from .run_all import pob
import os
import shutil
import unittest
from nose.tools.nontrivial import with_setup

class bestrefTest(unittest.TestCase):
    """ test for pob function in run_all.py
    """

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "pob_test_result")
        self.out_dir = os.path.join(self.test_dir, "plentyofbugs")
        self.data_dir = os.path.join(os.path.dirname(__file__), "test_data")
        self.plasmids_dir = os.path.join(self.data_dir, "plasmids", "")
        self.readsf = os.path.join(self.data_dir, "test_reads1.fq")
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
            
    def tearDown(self):
        "tear down test fixtures"
        shutil.rmtree(self.test_dir)


    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_pob(self):
        plasmids = (self.plasmids_dir)
        reads = (self.readsf)
        os.makedirs(self.test_dir)
        output_dir= (self.out_dir)
        test_result = pob(genomes_dir=plasmids, readsf=reads, output_dir=output_dir)
        print(test_result)
        assert "0.0340249" == test_result[1]
        assert os.path.basename("./test_data/plasmids/NC_009837.1.fasta") == \
            os.path.basename(test_result[0])
        return()

    


    
