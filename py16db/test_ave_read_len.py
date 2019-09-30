from .run_all import get_and_check_ave_read_len_from_fastq
import os
import shutil
import logging as logger
import unittest
import sys
import subprocess

class avereadlenTest(unittest.TestCase):
    ''' test for average read length
    '''
    def setUp(self):
        self.readsgunzipd = os.path.join(os.path.dirname(__file__), "test_data", "test_reads1.fq")
        self.readsgzipd = os.path.join(os.path.dirname(__file__), "test_data", "test_reads1.fq.gz")

    def tearDown(self):
        pass

    def test_get_ave(self):
        reads = self.readsgunzipd
        code, test_result = get_and_check_ave_read_len_from_fastq(fastq1=reads, logger=logger, minlen=50, maxlen=300)
        print(test_result)
        assert 150.0 == test_result
        assert code == 0

    def test_get_ave_too_long(self):
        reads = self.readsgunzipd
        code, test_result = get_and_check_ave_read_len_from_fastq(fastq1=reads, logger=logger, minlen=50, maxlen=100)
        print(test_result)
        assert 150.0 == test_result
        assert code == 2

    def test_get_ave_too_short(self):
        reads = self.readsgunzipd
        code, test_result = get_and_check_ave_read_len_from_fastq(fastq1=reads, logger=logger, minlen=200, maxlen=1000)
        print(test_result)
        assert 150.0 == test_result
        assert code == 1
