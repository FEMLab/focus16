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
# #16S from Escherichia coli genome {NC_011750.1.fna}

class test_alignmentData(unittest.TestCase):
    """ Using unittest for setUp and tearDown. This function creates the data
    necessary for alignment test"""

    def setUp(self):
        self.testdir = os.path.join(os.path.dirname(__file__), "test_data", "")
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
        
        
        self.shortbarrnap = os.path.join(os.path.dirname(__file__), "barrnapS")
        self.ecolibarrnap = os.path.join(os.path.dirname(__file__), "barrnapE")
        self.shortoutput = os.path.join(os.path.dirname(__file__), "test_data", "test_16s_multilineSHORT.fasta")
        self.ecoliout = os.path.join(os.path.dirname(__file__), "test_data", "ecoli")

        if os.path.exists(self.shortbarrnap):
            os.remove(self.shortbarrnap)
        if os.path.exists(self.ecolibarrnap):
            os.remove(self.ecolibarrnap)
        if os.path.exists(self.shortoutput):
            os.remove(self.shortoutput)
        if os.path.exists(self.ecoliout):
            os.remove(self.ecoliout)
        
    def tearDown(self):
        """ tear down test fixtures
        """
        if os.path.exists(self.shortbarrnap):
            os.remove(self.shortbarrnap)
        if os.path.exists(self.ecolibarrnap):
            os.remove(self.ecolibarrnap)
        if os.path.exists(self.shortbarrnap):
            os.remove(self.shortbarrnap)


    def test_alignment(self):
        accS = self.ecoliacc1
        ecolibarrnap = self.ecolibarrnap
        shortbarrnap = self.shortbarrnap
        shortoutput = self.shortoutput

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
           
        ## generates 16s sequences from 5 ecoli genomes
        ra.extract_16s_from_contigs(input_contigs=acc, barr_out=ecolibarrnap, output=ribo16, logger=logger)


        ## generates 16s sequences for 1 ecoli genome
        ra.extract_16s_from_contigs(input_contigs=accS, barr_out=shortbarrnap, output=shortoutput, logger=logger)

                       
        return()



## test_ave_read_len
## reads generated from {ART-bin-MountRainier-2016.06.05-MacOS64.tgz}
## Assuming art_bin_MountRainier/ is in ./16db/py16db 

class test_ave_read_length(unittest.TestCase):
    def setUp(self):
        """ Using the genome {NC_011750.1.fna}, produces reads for the 
        get_and_check_ave_read_len_from_fastq"""
        
        self.art = shutil.which("art_illumina")
        self.artreads = os.path.join(os.path.dirname(__file__), "test_data", "test_reads")
        self.genome = os.path.join(os.path.dirname(__file__), "test_data", "ecoli", "NC_011750.1.fna")
        
        if os.path.exists(self.artreads + "1.fq"):
            os.remove(self.artreads + "1.fq")
        if os.path.exists(self.artreads + "2.fq"):
            os.remove(self.artreads + "2.fq")


    def tearDown(self):
        alignmentfiles = glob.glob(os.path.join(self.artreads + "*.aln"))
        for file in alignmentfiles:
            os.remove(file)
            
                                                   
    def test_ave_read_len(self):
        artdir = self.art
        reads = self.artreads
        genome = self.genome = os.path.join(os.path.dirname(__file__), "test_data", "ecoli", "NC_011750.1.fna")

        cmd = "{artdir} -ss HS25 -i {genome} -p -l 150 -f 10 -m 400 -qs 10 -s 10 -o {reads}".format(**locals())

        subprocess.run(
            cmd,
            shell=sys.platform !="win32",
            check=True)
            
        return()

#sraFind file for test_get_sra_for_organism

class sraFind_Test(unittest.TestCase):
    #sraFind file from https://github.com/nickp60/sraFind/results/sraFind-All-biosample-with-SRA-hits.txt
    #or use fetch_sraFind_data.py
    
    def setUp(self):
        ''' generates mini sraFind file for test_get_sra_for_organism.py
        '''
        self.sraFind = os.path.join(os.path.dirname(__file__), "sraFind-All-biosample-with-SRA-hits.txt")
        self.test_sraFind = os.path.join(os.path.dirname(__file__), "test_data", "test_sraFind.txt")
        if os.path.exists(self.test_sraFind):
            os.remove(self.test_sraFind)

#taking 5 lines above and below target SRA for mini sraFind file
    def TestSmallsraFind(self):
        #grep command to take make small test_sraFind.txt
        sraFind = self.sraFind
        test_sraFind = self.test_sraFind
        org = "Lactobacillus oryzae" 
        cmd = "grep -C5 '{org}' {sraFind} > {test_sraFind}".format(**locals())
        
        subprocess.run(
            cmd,
            shell=sys.platform !="win32",
            check=True)


class parsefile_Test(unittest.TestCase):
    ''' writes the files necessary for the parse_status test
    '''
    def setUp(self):
        self.status = os.path.join(os.path.dirname(__file__), "test_data", "status", "")
        self.statusSRA = os.path.join(self.status, "sra")
        self.statusRIBO = os.path.join(self.status, "riboseed")

        if os.path.exists(self.status):
            shutil.rmtree(self.status)

    def test_writestatus(self):

        os.makedirs(self.status)

        with open(self.statusSRA, "a") as file1:
            file1.write("SRA COMPLETE\n")
        with open(self.statusRIBO, "a") as file2:
            file2.write("SRA COMPLETE\n")
            file2.write("RIBOSEED COMPLETE\n")
            file2.write("PROCESSED\n")
            
#list of SRAs             
class sralist_Test(unittest.TestCase):
    def setUp(self):
        self.sralist = os.path.join(os.path.dirname(__file__), "test_data", "test_sralist.txt")
        
        if os.path.exists(self.sralist):
            os.remove(self.sralist)

    def test_sralist(self):
        with open(self.sralist, "a") as file:
            file.write("ERX3310125\n")
            file.write("ERX3289350\n")
            file.write("ERX3289335\n")
            file.write("SRX2141371\n")



        

