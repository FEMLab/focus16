#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import shutil
import subprocess
import logging
import glob


class fasterqdumpError(Exception):
    pass


class FocusDBData(object):
    def __init__(self, dbdir=None, sraFind_data=None):
        self.dbdir = dbdir
        self.SRAs_manifest = os.path.join(self.dbdir, "SRAs_manifest.tab")
        self.sraFind_data = None
        self.SRAs = {}
        self.setup_if_needed()
        self.read_manifest()

    def setup_if_needed(self):
        """ if needed, create directory and write header for status file
        focusDB saves all data to a .focusDB dir in ones home folder
        """
        if not os.path.exists(self.dbdir):
            os.makedirs(self.dbdir)
        if not os.path.exists(self.SRAs_manifest):
            with open(self.SRAs_manifest, "w") as outf:
                outf.write("SRA_accession\tStatus\tOrganism\n")

    def read_manifest(self):
        with open(self.SRAs_manifest, "r") as inf:
            for i, line in enumerate(inf):
                acc, status, org = line.strip().split("\t")
                if i != 0:
                    self.SRAs[acc] = {
                        "status": status,
                        "organism": org,
                    }

    def update_manifest(self, newacc, newstatus, organism, logger):
        tmp = self.SRAs_manifest + ".bak"
        shutil.move(self.SRAs_manifest, tmp)
        with open(tmp, "r") as inf, open(self.SRAs_manifest, "w") as outf:
            for i, line in enumerate(inf):
                acc, status, lineorg = line.strip().split("\t")
                if acc == newacc:
                    pass
                else:
                    outf.write("{}\t{}\t{}\n".format(acc, status, lineorg))
            # if we still haven;t written out our new SRA (ie, if we are adding
            # a new one, not updating)
            outf.write("{}\t{}\t{}\n".format(newacc, newstatus, organism))
        os.remove(tmp)
        self.read_manifest()

    def fetch_sraFind_data(self, dest_path, logger):
        if dest_path is None:
            dest_path = os.path.join(
                self.dbdir, "sraFind-All-biosample-with-SRA-hits.txt")
        sraFind_results = str(
            "https://raw.githubusercontent.com/nickp60/sraFind/master/"+
            "results/sraFind-All-biosample-with-SRA-hits.txt"
        )
        # gets just the file name
        if not os.path.exists(dest_path):
            logger.info("Downloading sraFind Dump")
            download_sraFind_cmd = str("wget " + sraFind_results + " -O " + dest_path)
            logger.debug(download_sraFind_cmd)
            subprocess.run(
                download_sraFind_cmd,
                shell=sys.platform != "win32",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True)
        self.sraFind_data = dest_path

    def get_SRA_data(self, cores, SRA, org, logger):
        """download_SRA_if_needed
        This doesnt check the manifest right off the bad to make it easier for
        users to move data into the .focusdb dir manually
        1) check if dir and files exists.  If so, recheck and return files path
        if all is well;
        2) check manifest; there will be an status message if we removed the
        data intentionally.  raises an error if it looks like files have gone missing
        3) rerun if needed, and return the results
        """
        suboutput_dir_raw = os.path.join(self.dbdir, SRA, "")

        rawreadsf, rawreadsr, download_error_message = \
            self.check_fastq_dir(this_data=suboutput_dir_raw, logger=logger)
        if download_error_message == "":
            self.update_manifest(
                newacc=SRA,
                newstatus="PASS",
                organism=org,
                logger=logger)
            return (rawreadsf, rawreadsr, download_error_message)
        elif SRA in self.SRAs.keys():
            if self.SRAs[SRA]['status'] == "-":
                raise ValueError("Corrupted manifest; unable to find previously processed files")
            elif self.SRAs[SRA]['status'] == "DOWNLOAD ERROR":
                if os.path.exists(suboutput_dir_raw):
                    shutil.rmtree(suboutput_dir_raw)
            elif self.SRAs[SRA]['status'] == "LIBRARY TYPE ERROR":
                # dont try to reprocess
                return (None, None, "Library type Error")
        elif os.path.exists(suboutput_dir_raw):
                shutil.rmtree(suboutput_dir_raw)
        else:
            pass
        os.makedirs(suboutput_dir_raw)
        # defaults to 6 threads or whatever is convenient; we suspect I/O limits
        # using more in most cases, so we don't give the user the option
        # to increase this
        # https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
        cmd = str("fasterq-dump {SRA} -O " +
                  "{suboutput_dir_raw} --split-files").format(**locals())
        logger.info("Downloading %s", SRA)
        logger.debug(cmd)
        try:
            subprocess.run(cmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        except subprocess.CalledProcessError:
            self.update_manifest(
                newacc=SRA,
                newstatus="DOWNLOAD ERROR",
                organism=org,
                logger=logger)
            logger.critical("Error running fasterq-dump")
            raise fasterqdumpError
        rawreadsf, rawreadsr, download_error_message = \
            self.check_fastq_dir(this_data=suboutput_dir_raw, logger=logger)
        if download_error_message == "":
            self.update_manifest(
                newacc=SRA,
                newstatus="PASS",
                organism=org,
                logger=logger)
            return (rawreadsf, rawreadsr, download_error_message)
        else:
            return (None, None, download_error_message)

    def check_fastq_dir(self, this_data, logger):
        # Double check the download worked.  If its a single lib,
        # it wont have a _1 prefix, so we check if that exists
        # and if so, adjust expeactations
        if not os.path.exists(this_data):
            # this never happens unless a run is aborted;
            # regardless, we want to make sure we attempt to re-download
            return(None, None, "No directory created for SRA data during download")
        fastqs = glob.glob(os.path.join(this_data, "", "*.fastq"))
        logger.debug("fastqs detected: %s", " ".join(fastqs))
        if len(fastqs) == 0:
            return(None, None, "No fastq files downloaded")
        rawf, rawr, raws = [], [], []
        rawreadsf, rawreadsr = None, None
        download_error_message = ""
        for fastq in fastqs:
            if fastq.endswith("_1.fastq"):
                # not appending, cause we want to bump out any single libs t
                # that may have been read in first
                if len(rawf) != 0:
                    logger.warning("ignoring extra library %s", " ".join(rawf))
                rawf = [fastq]
            elif fastq.endswith("_2.fastq"):
                rawr.append(fastq)
            elif fastq.endswith(".fastq") and not fastq.endswith("_3.fastq"):
                if len(rawf) == 0:
                    # This is how we treat single libraries
                    rawf.append(fastq)
                else:
                    logger.warning("ignoring extra library %s", fastq)
            else:
                logger.error("Can only process paired or reads")
                logger.error(fastqs)
                download_error_message = "Unexpected item in the bagging area"
                download_error_message = "Unable to process library type"
        if len(set(rawf)) == 1:
            rawreadsf = rawf[0]
        elif len(set(rawf)) > 1:
            download_error_message = "multiple forward reads files detected"
        else:
            download_error_message = "No forward/single read files detected"

        if len(set(rawr)) == 1:
            rawreadsr = rawr[0]
        elif len(set(rawr)) > 1:
            download_error_message = "multiple reverse reads files detected"
        else:
            rawreadsr= None
        # catch only .fastq and _2.fastq weird combo
        if rawreadsf is not None:
             if not rawreadsf.endswith("_1.fastq") and rawreadsr is not None:
                 download_error_message = "cannot process a single library " + \
                     "file and a reverse file"
        return (rawreadsf, rawreadsr, download_error_message)
