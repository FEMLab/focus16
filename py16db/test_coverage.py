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

        self.reads1 = os.path.join(os.path.dirname(__file__), "test_data",  "test_reads1.fq")
        self.reads2 = os.path.join(os.path.dirname(__file__), "test_data",  "test_reads2.fq")
        
        self.downsample_dir = os.path.join(self.test_dir, "downsampled")
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def tearDown(self):
        "tear down test fixtures"
        for dir in [self.test_dir, self.downsample_dir]:
            if os.path.exists(dir):
                shutil.rmtree(dir)
    
    def test_coverage(self):
        #genome is from NC_011750.1 ~5132068bp at 10X coverage
        #reads are generated from this, under 16db/py16db/generator.py
        
        reads1 = self.reads1
        reads2 = self.reads2
        test_result = get_coverage(approx_length=5132068, fastq1=reads1, fastq2=reads2, read_length=150, logger=logger)
       
        assert 9.9997554592028 == test_result
        return()
    
    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_downsample_PE(self):
            #genome is built from NC_011750.1 ~5132068bp at 10X coverage
        #os.makedirs(self.downsample_dir, )
        reads1, reads2 = downsample(
            read_length=150,
            approx_length=5132068,
            fastq1=self.reads1,
            fastq2=self.reads2,
            destination=self.downsample_dir, 
            maxcoverage=2,
            logger=logger)
        down_cov1 = get_coverage(read_length=150, approx_length=5132068, fastq1=reads1, fastq2=reads2, logger=logger)
        down_cov2 = get_coverage(read_length=150, approx_length=5132068, fastq1=reads2, fastq2=reads2, logger=logger)
        assert 2.0110460344640795 == down_cov1
        assert 2.0110460344640795 == down_cov2
        return()
    
    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_downsample_SE(self):
        reads1, reads2 = downsample(
            read_length=150,
            approx_length=5132068,
            fastq1=self.reads1,
            fastq2=self.reads2,
            destination=self.downsample_dir, 
            maxcoverage=2,
            logger=logger)
        down_cov = get_coverage(read_length=150, approx_length=5132068, fastq1=reads1, fastq2=reads2, logger=logger)
        assert 2.0110460344640795 == down_cov
        return()
        
