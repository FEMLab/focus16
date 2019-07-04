from .run_all import alignment
from nose.tools.nontrivial import with_setup
import os
import shutil
import unittest
import logging as logger

class alignmentTest(unittest.TestCase):
   """ test for alignment function
   """
   def setUp(self):
      self.test_dir = os.path.join(os.path.dirname(__file__),
                                   "alignment_test")
      self.fasta = os.path.join(os.path.dirname(__file__),
                                "test_data", "test_16s_multiline.fasta")
      if os.path.exists(self.test_dir):
         shutil.rmtree(self.test_dir)
         
   def tearDown(self):
      """ tear down test fixtures
      """
      shutil.rmtree(self.test_dir)

   @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                    "skipping this test on travis.CI")    
   def test_alignment(self):
      test_output = (self.test_dir)
      os.makedirs(test_output)
      test_fasta = (self.fasta)
      test_result = alignment(fasta=test_fasta, output=test_output, logger=logger)
      print(test_result)
      mafftoutput = os.path.join(self.test_dir, "alignment", "MSA.fasta")
      with open(mafftoutput, "r") as infile:
         firstline = infile.readline().strip()
         print(firstline)
      assert firstline == ">CP003686.1 :16889:18441"
