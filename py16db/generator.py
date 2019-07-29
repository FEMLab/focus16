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

# ##test_alignment 
# #16S from Streptococcus mutans complete genome {CP003686.1}

class test_alignmentData(unittest.TestCase):
    """ Using unittest for setUp and tearDown. This function creates the data
    necessary for alignment test"""

    def setUp(self):
        self.testdir = os.path.join(os.path.dirname(__file__), "test_data", "")
        self.acc = os.path.join(self.testdir, "CP003686.1.fna")
        self.ecolidir = os.path.join(self.testdir, "ecoli")

        self.ecoliacc1 = os.path.join(self.ecolidir, "NC_011750.1.fna")
        self.ecoliurl1 = "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/026/345/GCA_000026345.1_ASM2634v1/GCA_000026345.1_ASM2634v1_genomic.fna.gz"

        self.ecoliacc2 = os.path.join(self.ecolidir, "NC_000913.3.fna")
        self.ecoliurl2 = "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/005/845/GCA_000005845.2_ASM584v2/GCA_000005845.2_ASM584v2_genomic.fna.gz"

        self.ecoliacc3 = os.path.join(self.ecolidir, "NC_017634.1.fna")
        self.ecoliurl3 = "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/183/345/GCA_000183345.1_ASM18334v1/GCA_000183345.1_ASM18334v1_genomic.fna.gz" 

        self.ecoliacc4 = os.path.join(self.ecolidir, "NC_018658.1.fna")
        self.ecoliurl4 = "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/299/455/GCA_000299455.1_ASM29945v1/GCA_000299455.1_ASM29945v1_genomic.fna.gz" 

        self.ecoliacc5 = os.path.join(self.ecolidir, "NC_011739.1.fna")
        self.ecoliurl5 = "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/026/325/GCA_000026325.2_ASM2632v2/GCA_000026325.2_ASM2632v2_genomic.fna.gz"
        
        
        self.strepbarrnap = os.path.join(os.path.dirname(__file__), "barrnapS")
        self.ecolibarrnap = os.path.join(os.path.dirname(__file__), "barrnapE")
        self.strepout = os.path.join(os.path.dirname(__file__), "test_data", "test_16s_multilineSHORT.fasta")
        self.ecoliout = os.path.join(os.path.dirname(__file__), "test_data", "ecoli")
        if os.path.exists(self.strepbarrnap):
            os.remove(self.strepbarrnap)
        if os.path.exists(self.ecolibarrnap):
            os.remove(self.ecolibarrnap)
        if os.path.exists(self.acc):
            os.remove(self.acc)
        if os.path.exists(self.strepout):
            os.remove(self.strepout)
        if os.path.exists(self.ecoliout):
            os.remove(self.ecoliout)
        
    def tearDown(self):
        """ tear down test fixtures
        """
        if os.path.exists(self.strepbarrnap):
            os.remove(self.strepbarrnap)
        if os.path.exists(self.ecolibarrnap):
            os.remove(self.ecolibarrnap)


    def test_alignment(self):
        acc = self.acc
        strepbarrnap = self.strepbarrnap
        ecolibarrnap = self.ecolibarrnap
        strepoutput = self.strepout

        WGETcmd =  "wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/271/865/GCA_000271865.1_ASM27186v1/GCA_000271865.1_ASM27186v1_genomic.fna.gz -O {acc}.gz".format(**locals())
        subprocess.run(
            WGETcmd,
            shell=sys.platform != "win32",
            check=True)
        
        GUNZIPcmd = "gunzip {acc}.gz".format(**locals())
        subprocess.run(
            GUNZIPcmd,
            shell=sys.platform != "win32",
            check=True)
        
        ra.extract_16s_from_contigs(input_contigs=acc, barr_out=strepbarrnap, output=strepoutput, logger=logger)

        ##Now for the 5 ecoli genomes
        os.makedirs(self.ecolidir)
        
        ecolicmd1 = "{self.ecoliurl1} -O {self.ecoliacc1}.gz".format(**locals())
        ecolicmd2 = "{self.ecoliurl2} -O {self.ecoliacc2}.gz".format(**locals())
        ecolicmd3 = "{self.ecoliurl3} -O {self.ecoliacc3}.gz".format(**locals())
        ecolicmd4 = "{self.ecoliurl4} -O {self.ecoliacc4}.gz".format(**locals())
        ecolicmd5 = "{self.ecoliurl5} -O {self.ecoliacc5}.gz".format(**locals())
        
        ribo16 = os.path.join(self.testdir, "ribo16")
        if os.path.exists(ribo16):
            os.remove(ribo16)
        
        for cmd in [ecolicmd1, ecolicmd2, ecolicmd3, ecolicmd4, ecolicmd5]:
            subprocess.run(cmd, 
                           shell=sys.platform != "win32",
                           check=True)

        for acc in [self.ecoliacc1, self.ecoliacc2, self.ecoliacc3, self.ecoliacc4, self.ecoliacc5]:
            gzipcmd = "gunzip {acc}.gz".format(**locals())
            subprocess.run(gzipcmd,
                           shell=sys.platform !="win32",
                           check=True)
           
        
        ra.extract_16s_from_contigs(input_contigs=acc, barr_out=ecolibarrnap, output=ribo16, logger=logger)
                           
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

        cmd = "{artdir} -ss HS25 -i {genome} -p -l 150 -f 10 -m 400 -s 10 -o {reads}".format(**locals())

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
     

