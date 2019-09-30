from .align_and_trim_focusdb import mafft as alignment
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
                                "test_data", "test_16s_multilineSHORT.fasta")

   def test_alignment(self):
      test_output = (self.test_dir)
      test_fasta = (self.fasta)
      outpath, test_result = alignment(multifasta=test_fasta, pre=test_output)
      assert test_result.startswith("mafft --retree 2 --reorder"), "malformatted mafft cmd"
      assert outpath == self.test_dir + ".mafft"


class alignmentLongTest(unittest.TestCase):
   """ test for the alignment step, using 5 whole genome sequences from E. coli
   """
   def setUp(self):
      self.test2_dir = os.path.join(os.path.dirname(__file__),
                                    "alignment_test_long")
      self.inputfasta = os.path.join(os.path.dirname(__file__), "test_data", "ribo16")

   def test_alignmentlong(self):
      test_output = (self.test2_dir)
      inputfasta = (self.inputfasta)
      out_path, test_result = alignment(multifasta=inputfasta, pre=test_output)
      assert test_result.startswith("mafft --retree 2 --reorder"), "malformatted mafft cmd"
      assert out_path == self.test2_dir + ".mafft"
