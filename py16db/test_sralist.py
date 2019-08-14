from .run_all import sralist
import unittest
import os

class sralistTest(unittest.TestCase):
    def setUp(self):
        self.sralist = os.path.join(os.path.dirname(__file__), "test_data", "test_sralist.txt")
        

    def test_sra_list(self):
        test_result = sralist(list=self.sralist)
        print(test_result)
        assert test_result == ["ERX3310125", "ERX3289350", "ERX3289335", "SRX2141371"]
