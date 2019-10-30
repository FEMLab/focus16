from .shared_methods import filter_sraFind, get_ave_read_len_from_fastq

from .run_focusDB import  check_read_len
import os
import shutil
import unittest
from nose.tools.nontrivial import with_setup
import logging as logger
import math

class filter_SRATest(unittest.TestCase):
    ''' test for filter_srapure and download_sra in run_all.py
    '''

    def setUp(self):
        self.sra_find=os.path.join(os.path.dirname(__file__), "test_data", "test_sraFind.txt")

    def test_filter_SRA(self):
        test_result = filter_sraFind(sraFind=self.sra_find,
                                 organism_name="Lactobacillus oryzae",
                                 thisseed=1,
                                 use_available=False,
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


class avereadlenTest(unittest.TestCase):
    ''' test for average read length
    '''
    def setUp(self):
        self.readsgunzipd = os.path.join(os.path.dirname(__file__), "test_data", "test_reads1.fq")
        self.readsgzipd = os.path.join(os.path.dirname(__file__), "test_data", "test_reads1.fq.gz")

    def test_get_ave(self):
        reads = self.readsgunzipd
        test_result = get_ave_read_len_from_fastq(
            fastq1=reads, logger=logger)
        assert 150 == math.floor(test_result)

    def test_get_ave_too_long(self):
        reads = self.readsgunzipd
        test_result = get_ave_read_len_from_fastq(
            fastq1=reads, logger=logger)
        code = check_read_len(test_result, logger=logger, minlen=50, maxlen=100)
        print(test_result)
        assert 150 == math.floor(test_result)
        assert code == 2

    def test_get_ave_too_short(self):
        reads = self.readsgunzipd
        test_result = get_ave_read_len_from_fastq(
            fastq1=reads, logger=logger)
        code = check_read_len(test_result, logger=logger, minlen=200, maxlen=1000)
        print(test_result)
        assert 150 == math.floor(test_result)
        assert code == 1
