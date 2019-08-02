from .run_all import run_riboseed
import os
import shutil
import unittest
import logging as logger
from nose.tools.nontrivial import with_setup


class coverageTests(unittest.TestCase):
    """ tests for coverage and downsample functions in run_all.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "riboSeed")
        self.readsf = os.path.join(os.path.dirname(__file__), "test_data",
                                   "reads1.fq")
        self.readsr = os.path.join(os.path.dirname(__file__), "test_data",
                                   "reads2.fq")
        self.sra = os.path.join(os.path.dirname(__file__), "test_data",
                                "ecoli", "NC_011750.1.fna")
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def tearDown(self):
        "tear down test fixtures"
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
            
    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")    
    def test_riboseed(self):
        readsf = (self.readsf)
        readsr = (self.readsr)
        output_dir = (self.test_dir)
        os.makedir = (output_dir)
        sra = (self.sra)
        
        test_result = run_riboseed(sra=sra, readsf=readsf,
                                   readsr=readsr, cores="4", threads="1",
                                   subassembler="spades",
                                   memory=8,
                                   output=output_dir, logger=logger)

        assert test_result == "ribo run -r /Users/alexandranolan/Desktop/16db/py16db/test_data/ecoli/NC_011750.1.fna -F /Users/alexandranolan/Desktop/16db/py16db/test_data/reads1.fq -R /Users/alexandranolan/Desktop/16db/py16db/test_data/reads2.fq --cores 4 --threads 1 -v 1 --serialize -o /Users/alexandranolan/Desktop/16db/py16db/riboSeed --subassembler spades --stages score --memory 8"
    
