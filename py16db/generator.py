from . import run_all as ra
import glob
import os
import sys
import unittest
import shutil
import subprocess
import logging as logger
import argparse

class test_requirements(unittest.TestCase):
    """ Ensuring the correct tools are present
    """
    def test_check_requirements(self):
        if not shutil.which("art_illumina"):
            logger.debug("art_bin_MountRainier not cloned, see README.md")
            sys.exit(1)

##test_alignment 
#16S from Streptococcus mutans complete genome {CP003686.1}

class test_alignmentData(unittest.TestCase):
    """ Using unittest for setUp and tearDown. This function creates the data
    necessary for alignment test"""

    def setUp(self):
        self.genome = os.path.join(os.path.dirname(__file__), "test_data", "CP003686.1.fna")
        self.barrnap = os.path.join(os.path.dirname(__file__), "barrnap")
        self.out = os.path.join(os.path.dirname(__file__), "test_data", "test_16s_multiline.fasta")
        if os.path.exists(self.barrnap):
            os.remove(self.barrnap)
        if os.path.exists(self.genome):
            os.remove(self.genome)
         
    def tearDown(self):
        """ tear down test fixtures
        """
        if os.path.exists(self.barrnap):
            os.remove(self.barrnap)

    def test_alignment(self):
        genome = self.genome
        barrnap = self.barrnap
        output = self.out
       

        WGETcmd =  "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/271/865/GCA_000271865.1_ASM27186v1/GCA_000271865.1_ASM27186v1_genomic.fna.gz -O {genome}.gz".format(**locals())
        subprocess.run(
            WGETcmd,
            shell=sys.platform != "win32",
            check=True)
        
        GUNZIPcmd = "gunzip {genome}.gz".format(**locals())
        subprocess.run(
            GUNZIPcmd,
            shell=sys.platform != "win32",
            check=True)
        
        ra.extract_16s_from_contigs(input_contigs=genome, barr_out=barrnap, output=output, logger=logger)
        return()

## test_ave_read_len
## reads generated from {ART-bin-MountRainier-2016.06.05-MacOS64.tgz}
## Assuming art_bin_MountRainier/ is in ./16db/py16db 

class test_ave_read_length(unittest.TestCase):
    def setUp(self):
        """ Using the genome {CP003686}, produces reads for the 
        get_and_check_ave_read_len_from_fastq"""
        self.art = shutil.which("art_illumina")
        self.artreads = os.path.join(os.path.dirname(__file__), "test_data", "test_reads")
        self.genome = os.path.join(os.path.dirname(__file__), "test_data", "CP003686.1.fna")
        
        if os.path.exists(self.artreads):
            os.remove(self.artreads)

    def tearDown(self):
        alignmentfiles = glob.glob(os.path.join(self.artreads + "*.aln"))
        for file in alignmentfiles:
            os.remove(file)
            
        
                                                     
    def test_ave_read_len(self):
        artdir = self.art
        reads = self.artreads
        genome = self.genome = os.path.join(os.path.dirname(__file__), "test_data", "CP003686.1.fna")

        cmd = "{artdir} -ss HS25 -i {genome} -p -l 150 -f 20 -m 400 -s 10 -o {reads}".format(**locals())

        subprocess.run(
            cmd,
            shell=sys.platform !="win32",
            check=True)
            
        return()

##test_best_ref
##plasmids are obtained from 

# class test_bestRef(unittest.TestCase):
#     def setUp(self):
#         self.reads = self.artreads = os.path.join(os.path.dirname(__file__), "test_data", "test_reads1.fq")
#         self.plasmids = 
     

