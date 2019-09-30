from .run_all import make_riboseed_cmd
import os
import shutil
import unittest
import sys
import subprocess
import logging as logger
from nose.tools.nontrivial import with_setup


class coverageTests(unittest.TestCase):
    """ tests for coverage and downsample functions in run_all.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "riboSeed")
        self.data_dir = os.path.join(os.path.dirname(__file__), "test_data", "")
        self.readsgunzipd1 = os.path.join(self.data_dir, "test_reads1.fq")
#        self.readsgzipd1 = os.path.join(self.data_dir, "test_reads1.fq.gz")

        self.readsgunzipd2 = os.path.join(self.data_dir, "test_reads2.fq")
#        self.readsgzipd2 = os.path.join(self.data_dir, "test_reads2.fq.gz")
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
        readsf = self.readsgunzipd1
        readsr = self.readsgunzipd2
        output_dir = self.test_dir
        os.makedir = output_dir
        sra = (self.sra)

        test_result = make_riboseed_cmd(sra=sra, readsf=readsf,
                                   readsr=readsr, cores="4", threads="1",
                                   subassembler="spades",
                                   memory=8,
                                   output=output_dir, logger=logger)
        target_cmd = "ribo run -r /Users/alexandranolan/Desktop/16db/py16db/test_data/ecoli/NC_011750.1.fna -F /Users/alexandranolan/Desktop/16db/py16db/test_data/test_reads1.fq -R /Users/alexandranolan/Desktop/16db/py16db/test_data/test_reads2.fq --cores 4 --threads 1 -v 1 -o /Users/alexandranolan/Desktop/16db/py16db/riboSeed --serialize --subassembler spades --just_seed --memory 8"
        for part in range(len(target_cmd.split(" "))):
            if part not in [3, 5, 7, 15]:
                print(test_result.split(" ")[part] )
                print(target_cmd.split(" ")[part] )
                assert test_result.split(" ")[part] == target_cmd.split(" ")[part]
