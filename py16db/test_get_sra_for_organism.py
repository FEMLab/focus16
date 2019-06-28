from .run_all import download_SRA, filter_srapure
import os
import shutil
import unittest
from nose.tools.nontrivial import with_setup

class get_sra_for_organismTest(unittest.TestCase):
    ''' test for filter_srapure and download_sra in run_all.py'''
    def setUp(self):
        self.test_dir=os.path.join(os.path.dirname(__file__), "test_function", "")
        self.sra_find=os.path.join(os.path.dirname(__file__), "test_data", "test_sraFind.txt")
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
            
    def teardown(self):
        "tear down test fixtures"
        shutil.rmtree(self.test_dir)
        
    def test_filter_srapure(self):
        test_result = filter_srapure(path=self.sra_find,
                                     organism_name="Lactobacillus oryzae", 
                                     strains=1, get_all=True)
        assert ["DRR021662"] == test_result
            
    def test_download_SRA(self):
        os.makedirs(self.test_dir)
        download_SRA(destination=self.test_dir,
                     SRA="SRR8443698")
        
        assert os.path.isdir(self.test_dir), "not creating correct subdirectory"
                
