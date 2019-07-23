from .run_all import get_coverage, downsample
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
                                     "downsample_test_result")

        self.reads = os.path.join(os.path.dirname(__file__), "test_data",  "test_reads1.fq")
        self.downsample_dir = os.path.join(self.test_dir, "downsampled")
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def tearDown(self):
        "tear down test fixtures"
        for dir in [self.test_dir, self.downsample_dir]:
            if os.path.exists(dir):
                shutil.rmtree(dir)
    
    def test_coverage(self):
        #genome is built from NC_0176759.1 ~150000bp
        
        reads = (self.reads)
        test_result = get_coverage(approx_length=150000, fastq1=reads, read_length=150, logger=logger)
        print(test_result)
        assert 9.79975  == test_result
        return()
    
    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_downsample_PE(self):
            #genome is built from NC_0176759.1 ~150000bp
        os.makedirs(self.downsample_dir)
        reads1, reads2 = downsample(
            read_length=150,
            approx_length=150000,
            fastq1=self.reads,
            fastq2=self.reads,
            destination=self.test_dir, 
            maxcoverage=2,
            logger=logger)
        down_cov1 = get_coverage(read_length=150, approx_length=150000, fastq1=reads1, logger=logger)
        down_cov2 = get_coverage(read_length=150, approx_length=150000, fastq1=reads2, logger=logger)
        assert 1.99275 == down_cov1
        assert 1.99275 == down_cov2
        return()
    
    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_downsample_SE(self):
        os.makedirs(self.downsample_dir)
        reads1, reads2 = downsample(
            read_length=150,
            approx_length=150000,
            fastq1=self.reads,
            fastq2=None,
            destination=self.test_dir, 
            maxcoverage=2,
            logger=logger)
        down_cov = get_coverage(read_length=150, approx_length=150000, fastq1=reads1, logger=logger)
        assert 1.99275 == down_cov
        return()
        
