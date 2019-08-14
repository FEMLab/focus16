from .run_all import parse_status_file
from .run_all import main
import os 
import unittest
import shutil
from nose.tools.nontrivial import with_setup



class parse_statusTest(unittest.TestCase):
    def setUp(self):
        ##sra is a status file with just sra checkpoint
        ##riboseed is a status file with just riboseed checkpoint    

        self.statusSRA = os.path.join(os.path.dirname(__file__), "test_data",  "status", "sra")
        self.statusriboseed = os.path.join(os.path.dirname(__file__), "test_data", "status", "riboseed")

    
    def test_parse_status_SRA(self):
        statusfile = (self.statusSRA)
        test_result = parse_status_file(path = statusfile)
        print(test_result)
        assert ["SRA COMPLETE"] == test_result

    def test_parse_status_RIBO(self):
        statusfile = (self.statusriboseed)
        test_result = parse_status_file(path = statusfile)
        assert ["SRA COMPLETE", "RIBOSEED COMPLETE", "PROCESSED"] == test_result

    
        
        
    
