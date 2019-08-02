from .run_all import check_rDNA_copy_number
import os
import unittest
import shutil
import logging as logger

class checkrDNATest(unittest.TestCase):
    def setUp(self):
        self.test_ref = os.path.join(os.path.dirname(__file__), "test_data", "ecoli", "NC_011750.1.fna")
        self.test_dir = os.path.join(os.path.dirname(__file__), "check_test")
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def tearDown(self):
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_check_rDNA(self):
        # NC_011750.1.fna has 7 rrn operons
        ref = (self.test_ref)
        out = (self.test_dir)
        os.makedirs(out)
        test_result = check_rDNA_copy_number(ref=ref, output=out, logger=logger)
        assert test_result == 7
        
