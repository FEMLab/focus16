from .run_all import filter_SRA
import os
import shutil
import unittest
from nose.tools.nontrivial import with_setup
import logging as logger

class filter_SRATest(unittest.TestCase):
    ''' test for filter_srapure and download_sra in run_all.py
    '''

    def setUp(self):
        self.sra_find=os.path.join(os.path.dirname(__file__), "test_data", "test_sraFind.txt")

    def test_filter_SRA(self):
        test_result = filter_SRA(sraFind=self.sra_find,
                                     organism_name="Lactobacillus oryzae",
                                     thisseed=1,
                                     strains=1, get_all=True, logger=logger)
        assert ["DRR021662"] == test_result


class download_SRATest(unittest.TestCase):
    ''' test for filter_srapure and download_sra in run_all.py
    '''

    def setUp(self):
        self.test_dir=os.path.join(os.path.dirname(__file__), "test_function", "")
        self.sra_find=os.path.join(os.path.dirname(__file__), "test_data", "test_sraFind.txt")
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)

    def tearDown(self):
        "tear down test fixtures"
        shutil.rmtree(self.test_dir)

    # @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
    #                  "skipping this test on travis.CI")
    # def test_download_SRA(self):

    #     download_SRA(cores=4, destination=self.test_dir,
    #                  SRA="SRR8443698", logger=logger)

    #     assert os.path.exists(os.path.join(self.test_dir, "SRR8443698_1.fastq"))
