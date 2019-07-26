from .run_all import extract_16s_from_contigs
from nose.tools.nontrivial import with_setup
import os
import shutil
import unittest
import logging as logger

class extractTest(unittest.TestCase):
   """ test for extract_16s_from_contigs 
   """
   def setUp(self):
      self.test_dir = os.path.join(os.path.dirname(__file__),                     
                                   "extract_test")
      self.contigs = os.path.join(os.path.dirname(__file__),
                                  "test_data", "ecoli", "NC_000913.3.fna")
      self.barrnap = os.path.join(self.test_dir, "barrnap")
      self.out_dir = os.path.join(self.test_dir, "ribo16")
      if os.path.exists(self.test_dir):
         shutil.rmtree(self.test_dir)
      
   def tearDown(self):
      """tear down test fixtures
      """
      shutil.rmtree(self.test_dir)

   @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                    "skipping this test on travis.CI")
   def test_extract(self):
      os.makedirs(self.test_dir)
      input_contigs = (self.contigs)
      barr_out = (self.barrnap)
      output = (self.out_dir)
      test_result=extract_16s_from_contigs(input_contigs=input_contigs,
                                           barr_out=barr_out, output=output,
                                           logger=logger)
      assert os.path.exists(test_result)
       

class extractTestLong(unittest.TestCase):
   """ test for extract_16s_from_contigs 
   """
   def setUp(self):
      self.test_dir = os.path.join(os.path.dirname(__file__),                     
                                   "extract_test")
      self.contigs = os.path.join(os.path.dirname(__file__),
                                  "test_data", "ecoli", "")
      self.barrnap = os.path.join(self.test_dir, "barrnap")
      self.out_dir = os.path.join(self.test_dir, "extracttestribo16")
      if os.path.exists(self.test_dir):
         shutil.rmtree(self.test_dir)
      
   def tearDown(self):
      """tear down test fixtures
      """
      shutil.rmtree(self.test_dir)

   @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                    "skipping this test on travis.CI")
   
   def test_extractLong(self):
      os.makedirs(self.test_dir)
      barr_out = (self.barrnap)
      output = (self.out_dir)

      for file in os.listdir(self.contigs):
         contig = os.path.join(self.contigs, file)
         test_result = extract_16s_from_contigs(input_contigs=contig,
                                                barr_out=barr_out, output=output,
                                                logger=logger)
         assert os.path.exists(output)


