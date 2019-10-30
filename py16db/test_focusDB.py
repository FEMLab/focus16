from .run_focusDB import get_coverage, downsample, make_riboseed_cmd, sralist,\
    pob, referenceNotGoodEnoughError
import os
import shutil
import unittest
import subprocess
import sys
import logging as logger
from nose.tools.nontrivial import with_setup


class coverageTests(unittest.TestCase):
    """ tests for coverage and downsample functions in run_all.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "downsample_test_result")

        self.data_dir = os.path.join(os.path.dirname(__file__), "test_data", "")

        self.readsgunzipd1 = os.path.join(self.data_dir, "test_reads1.fq")
        self.readsgzipd1 = os.path.join(self.data_dir, "test_reads1.fq.gz")

        self.readsgunzipd2 = os.path.join(self.data_dir, "test_reads1.fq")
        self.readsgzipd2 = os.path.join(self.data_dir, "test_reads1.fq.gz")
        self.downsample_dir = os.path.join(self.test_dir, "downsampled")


    def tearDown(self):
        "tear down test fixtures"
        for dir in [self.test_dir, self.downsample_dir]:
            if os.path.exists(dir):
                shutil.rmtree(dir)


    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_coverage(self):
        #genome is from NC_011750.1 ~5132068bp at 10X coverage
        #reads are generated from this, under 16db/py16db/generator.py

        reads1 = self.readsgunzipd1
        reads2 = self.readsgunzipd2
        test_result = get_coverage(approx_length=5132068, fastq1=reads1, fastq2=reads2, read_length=150, logger=logger)
        print(test_result)
        assert round(1.00, 2) == round(test_result, 2)
        return()

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_downsample_PE(self):
            #genome is built from NC_011750.1 ~5132068bp at 10X coverage
        #os.makedirs(self.downsample_dir, )
        reads1, reads2 = downsample(
            read_length=150,
            approx_length=5132068,
            fastq1=self.readsgunzipd1,
            fastq2=self.readsgunzipd2,
            destination=self.downsample_dir,
            mincoverage=.5,
            maxcoverage=2,
            run=True,
            logger=logger)
        down_cov1 = get_coverage(read_length=150, approx_length=5132068, fastq1=reads1, fastq2=reads2, logger=logger)
        print(down_cov1)
        # down_cov2 = get_coverage(read_length=150, approx_length=5132068, fastq1=reads2, fastq2=reads2, logger=logger)
        assert round(1.0110460344640795, 1) == round(down_cov1, 1)

        # assert 2.0110460344640795 == down_cov2

    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_downsample_SE(self):
        reads1, reads2 = downsample(
            read_length=150,
            approx_length=5132068,
            fastq1=self.readsgunzipd1,
            fastq2=self.readsgunzipd2,
            destination=self.downsample_dir,
            mincoverage=.5,
            maxcoverage=2,
            run=True,
            logger=logger)
        down_cov = get_coverage(read_length=150, approx_length=5132068, fastq1=reads1, fastq2=reads2, logger=logger)
        print(down_cov)
        assert round(1.00, 2 )  == round(down_cov, 2)


class sralistTest(unittest.TestCase):
    def setUp(self):
        self.sralist = os.path.join(os.path.dirname(__file__), "test_data", "test_sralist.txt")


    def test_sra_list(self):
        test_result = sralist(list=self.sralist)
        print(test_result)
        assert test_result == ["ERX3310125", "ERX3289350", "ERX3289335", "SRX2141371"]


class coverageTests(unittest.TestCase):
    """ tests for coverage and downsample functions in run_all.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "riboSeed")
        self.data_dir = os.path.join(os.path.dirname(__file__), "test_data", "")
        self.readsgunzipd1 = os.path.join(self.data_dir, "test_reads1.fq")
#        self.readsgzipd1 = os.path.join(self.data_dir, "test_reads1.fq.gz")

        self.readsgunzipd2 = os.path.join(self.data_dir, "test_reads2.fq")
#        self.readsgzipd2 = os.path.join(self.data_dir, "test_reads2.fq.gz")
        self.sra = os.path.join(os.path.dirname(__file__), "test_data",
                                "ecoli", "NC_011750.1.fna")
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)


    def tearDown(self):
        "tear down test fixtures"
        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)


    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_riboseed(self):
        readsf = self.readsgunzipd1
        readsr = self.readsgunzipd2
        output_dir = self.test_dir
        os.makedir = output_dir
        sra = (self.sra)

        test_result = make_riboseed_cmd(sra=sra, readsf=readsf,
                                        readsr=readsr, cores="4", threads="1",
                                        subassembler="spades",
                                        memory=8, just_seed=True,
                                        sge=False,
                                        output=output_dir, logger=logger)
        target_cmd = "ribo run -r /Users/alexandranolan/Desktop/16db/py16db/test_data/ecoli/NC_011750.1.fna -F /Users/alexandranolan/Desktop/16db/py16db/test_data/test_reads1.fq -R /Users/alexandranolan/Desktop/16db/py16db/test_data/test_reads2.fq --cores 4 --threads 1 -v 1 -o /Users/alexandranolan/Desktop/16db/py16db/riboSeed --serialize --subassembler spades --just_seed --skip_control --stages none --memory 8"
        for part in range(len(target_cmd.split(" "))):
            if part not in [3, 5, 7, 15]:
                print(test_result.split(" ")[part] )
                print(target_cmd.split(" ")[part] )
                assert test_result.split(" ")[part] == target_cmd.split(" ")[part]



class bestrefTest(unittest.TestCase):
    """ test for pob function in run_all.py
    """
    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(__file__),
                                     "pob_test_result", "")
        self.out_dir = os.path.join(self.test_dir, "plentyofbugs")
        self.data_dir = os.path.join(os.path.dirname(__file__), "test_data")
        self.plasmids_dir = os.path.join(self.data_dir, "ecoli", "")
        self.readsgunzipd = os.path.join(self.data_dir, "test_reads1.fq")
        self.readsgzipd = os.path.join(self.data_dir, "test_reads1.fq.gz")

        if os.path.exists(self.test_dir):
            shutil.rmtree(self.test_dir)
        if os.path.exists(os.path.join(self.plasmids_dir, "reference.msh")):
            os.remove(os.path.join(self.plasmids_dir, "reference.msh"))


    def tearDown(self):
        "tear down test fixtures"
        shutil.rmtree(self.test_dir)


    @unittest.skipIf("TRAVIS" in os.environ and os.environ["TRAVIS"] == "true",
                     "skipping this test on travis.CI")
    def test_pob(self):
        plasmids = (self.plasmids_dir)
        reads = (self.readsgunzipd)
        os.makedirs(self.test_dir)
        output_dir= (self.out_dir)
        with self.assertRaises(referenceNotGoodEnoughError):
            bad_test_result = pob(genomes_dir=plasmids, readsf=reads, output_dir=output_dir, maxdist=.05, logger=logger)
        test_result = pob(genomes_dir=plasmids, readsf=reads, output_dir=output_dir + "2", maxdist=.3, logger=logger)
        print(test_result)
        assert round(0.295981, 2) == round(test_result[1], 2)
