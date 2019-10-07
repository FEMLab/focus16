#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import shutil
import subprocess
import glob

from plentyofbugs import get_n_genomes as gng


class fasterqdumpError(Exception):
    pass


class FocusDBData(object):
    def __init__(self, dbdir=None, refdir=None,
                 sraFind_data=None, prokaryotes=None):
        self.dbdir = dbdir
        self.refdir = refdir
        self.prokaryotes = prokaryotes
        self.sraFind_data = sraFind_data
        self.SRAs = {}
        # get/set location of data
        self.get_focusDB_dir()
        # make dirs/files as needed
        self.refs_manifest = os.path.join(self.dbdir, "reference_genomes.tab")
        self.SRAs_manifest = os.path.join(self.dbdir, "SRAs_manifest.tab")
        self.setup_if_needed()
        self.read_SRA_manifest()

    def setup_if_needed(self):
        """ if needed, create directory and write header for status file
        focusDB saves all data to a .focusDB dir in ones home folder
        """
        if not os.path.exists(self.dbdir):
            os.makedirs(self.dbdir)
        if not os.path.exists(self.SRAs_manifest):
            with open(self.SRAs_manifest, "w") as outf:
                outf.write("SRA_accession\tStatus\tOrganism\n")

    def get_focusDB_dir(self):
        if self.dbdir is None:
            self.dbdir = os.path.join(os.path.expanduser("~"), ".focusDB", "")

    def check_genomes_dir(self, org):
        if org is None:
            raise ValueError("organism is required")
        if self.refdir is None:
            dirname = org.replace(" ", "_")
            self.refdir = os.path.join(
                self.dbdir, "references", dirname, "")
        else:
            # make sure we have a trailing pathsep for globs down the line
            self.refdir = os.path.join(self.refdir, "")
        if not os.path.exists(self.refdir):
            os.makedirs(self.refdir)

    def read_SRA_manifest(self):
        with open(self.SRAs_manifest, "r") as inf:
            for i, line in enumerate(inf):
                try:
                    acc, status, org = line.strip().split("\t")
                except ValueError as e:
                    print(self.SRAs_manifest)
                    print(line)
                    raise e
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
            # if we still haven't written out our new SRA (ie, if we are adding
            # a new one, not updating)
            outf.write("{}\t{}\t{}\n".format(newacc, newstatus, organism))
        os.remove(tmp)
        self.read_SRA_manifest()

    def fetch_sraFind_data(self, logger):
        if self.sraFind_data is None:
            self.sraFind_data = os.path.join(
                self.dbdir, "sraFind-All-biosample-with-SRA-hits.txt")
        sraFind_results = str(
            "https://raw.githubusercontent.com/nickp60/sraFind/master/" +
            "results/sraFind-All-biosample-with-SRA-hits.txt"
        )
        # gets just the file name
        if not os.path.exists(self.sraFind_data):
            logger.info("Downloading sraFind Dump")
            download_sraFind_cmd = str(
                "wget " + sraFind_results + " -O " + self.sraFind_data)
            logger.debug(download_sraFind_cmd)
            subprocess.run(
                download_sraFind_cmd,
                shell=sys.platform != "win32",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True)

    def run_prefetch_data(self, SRA_list, org, logger):

        pass

    def get_SRA_data(self, SRA, org, logger):
        """download_SRA_if_needed
        This doesnt check the manifest right off the bad to make it easier for
        users to move data into the .focusdb dir manually
        1) check if dir and files exists.  If so, recheck and return files path
        if all is well;
        2) check manifest; there will be an status message if we removed the
        data intentionally.
        raises an error if it looks like files have gone missing
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
                raise ValueError("Corrupted manifest; unable to find " +
                                 "previously processed files")
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
        # defaults to 6 threads or whatever is convenient;
        # we suspect I/O limits using more in most cases,
        # so we don't give the user the option to increase this
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
            self.update_manifest(
                newacc=SRA,
                newstatus="LIBRARY TYPE ERROR",
                organism=org,
                logger=logger)
            return (None, None, download_error_message)

    def check_fastq_dir(self, this_data, logger):
        # Double check the download worked.  If its a single lib,
        # it wont have a _1 prefix, so we check if that exists
        # and if so, adjust expeactations
        if not os.path.exists(this_data):
            # this never happens unless a run is aborted;
            # regardless, we want to make sure we attempt to re-download
            return(None, None, "No directory created for SRA during download")
        fastqs = glob.glob(os.path.join(this_data, "", "*.fastq"))
        logger.debug("fastqs detected: %s", " ".join(fastqs))
        if len(fastqs) == 0:
            return(None, None, "No fastq files downloaded")
        rawf, rawr = [], []
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
                logger.error("Can only process paired or single-end reads")
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
            rawreadsr = None
        # catch only .fastq and _2.fastq weird combo
        if rawreadsf is not None:
            if not rawreadsf.endswith("_1.fastq") and rawreadsr is not None:
                download_error_message = "cannot process a single library " + \
                    "file and a reverse file"
        return (rawreadsf, rawreadsr, download_error_message)

    #####################   Methods for dealing with reference genomes ########
    def decide_skip_or_download_genomes(self, args, logger):
        """ checks the genomes dir
        returns
        - 0 if all is well
        - 1 if no matching genomes in prokaryotes.txt
        - 2 if error downloading
        - 3 if error unzipping
        """
        # check basic integrity of genomes dir
        # it might exist, but making it here simplifies the control flow later
        # it seems counter intuitive, but checking the dir we might have just
        # created is easier than checking if it exists/is intact twice
        os.makedirs(self.refdir, exist_ok=True)
        if len(glob.glob(os.path.join(self.refdir, "*.fna"))) == 0:
            logger.info('Downloading genomes')
            return(self.our_get_n_genomes(
                org=args.organism_name,
                nstrains=args.nstrains,
                thisseed=args.seed,
                logger=logger)
            )
        if len(glob.glob(os.path.join(self.refdir, "*.fna.gz"))) != 0:
            logger.warning('Genome downloading may have been interupted; ' +
                           'downloading fresh')
            shutil.rmtree(self.refdir)
            os.makedirs(self.refdir)
            return(self.our_get_n_genomes(
                org=args.organism_name,
                nstrains=args.nstrains,
                thisseed=args.seed,
                logger=logger)
            )
        return 0

    def our_get_n_genomes(self, org, nstrains,  thisseed, logger):
        # taken from the main function of get_n_genomes
        if self.prokaryotes is None:
            self.prokaryotes = os.path.join(self.dbdir, "prokaryotes.txt")
        if not os.path.exists(self.prokaryotes):
            gng.fetch_prokaryotes(dest=self.prokaryotes)
        org_lines = gng.get_lines_of_interest_from_proks(path=self.prokaryotes,
                                                         org=org)
        if len(org_lines) == 0:
            return 1
        if nstrains == 0:
            nstrains = len(org_lines)
        print(nstrains)
        cmds = gng.make_fetch_cmds(
            org_lines,
            nstrains=nstrains,
            thisseed=thisseed,
            genomes_dir=self.refdir,
            SHUFFLE=True)
        print(cmds)
        # this can happen if tabs end up in metadata (see
        # AP019724.1 B. unifomis, and others
        # I updated plentyofbugs to try to catch it, but still might happen
        if len(cmds) == 0:
            return 1
        try:
            for i, cmd in enumerate(cmds):
                sys.stderr.write("Downloading genome %i of %i\n%s\n" %
                                 (i + 1, len(cmds), cmd))
                logger.debug(cmd)
                subprocess.run(
                    cmd,
                    shell=sys.platform != "win32",
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=True)
        except Exception as e:
            logger.error(e)
            return 2

        unzip_cmd = "gunzip " + os.path.join(self.refdir, "*.gz")
        sys.stderr.write(unzip_cmd + "\n")
        try:
            subprocess.run(
                unzip_cmd,
                shell=sys.platform != "win32",
                check=True)
        except Exception as e:
            logger.error(e)
            return 3
        return 0
