from .run_all import pob
import os
import shutil
import unittest
from nose.tools.nontrivial import with_setup
import logging as logger

class bestrefTest(unittest.TestCase):
    """ test for pob function in run_all.py
    """

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "pob_test_result")
        self.out_dir = os.path.join(self.test_dir, "plentyofbugs")
        self.data_dir = os.path.join(os.path.dirname(__file__), "test_data")
        self.plasmids_dir = os.path.join(self.data_dir, "ecoli", "")
        self.readsf = os.path.join(self.data_dir, "test_reads1.fq")
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
        if os.path.exists(os.path.join(self.plasmids_dir, "reference.msh")):
            os.remove(os.path.join(self.plasmids_dir, "reference.msh"))


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
        test_result = pob(genomes_dir=plasmids, readsf=reads, output_dir=output_dir, logger=logger)
        print(test_result)
        assert "0.00621379" == test_result[1]
        assert os.path.basename(os.path.join(self.plasmids_dir, "NC_011750.1.fna")) == \
            os.path.basename(test_result[0])
        return()

    


    
