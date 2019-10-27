#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import sys
import time
import shutil
import subprocess
import glob
import sqlite3
from plentyofbugs import get_n_genomes as gng


class fasterqdumpError(Exception):
    pass


class FocusDBData(object):
    def __init__(self, dbdir=None, refdir=None,
                 sraFind_data=None, prokaryotes=None,
                 krakendir=None):
        self.dbdir = dbdir
        self.refdir = refdir
        self.prokaryotes = prokaryotes
        self.sraFind_data = sraFind_data
        self.SRAs = {}  # acc, status, genus species
        self.krakendir = krakendir
        # get/set location of data
        self.get_focusDB_dir()
        # make dirs/files as needed
        self.refs_manifest = os.path.join(self.dbdir, "reference_genomes.tab")
        self.SRAs_manifest = os.path.join(self.dbdir, "SRAs_manifest.tab")
        self.db = os.path.join(self.dbdir, "focus.db")
        self.setup_if_needed()
        self.read_SRA_manifest()

    def setup_if_needed(self):
        """ if needed, create directory and write header for status file
        focusDB saves all data to a .focusDB dir in ones home folder
        """
        if not os.path.exists(self.dbdir):
            os.makedirs(self.dbdir)
        conn = sqlite3.connect(self.db)
        c = conn.cursor()
        # TODO store read len in DB
        # c.execute("CREATE TABLE IF NOT EXISTS SRAs (accession text PRIMARY KEY, status text, genus text, species text, readlen integer )")
        c.execute("CREATE TABLE IF NOT EXISTS SRAs (accession text PRIMARY KEY, status text, genus text, species text )")
        c.execute("CREATE TABLE IF NOT EXISTS Genomes (accssions text PRIMARY KEY, status text, genus text, species text)")
        # if not os.path.exists(self.SRAs_manifest):
        #     with open(self.SRAs_manifest, "w") as outf:
        #         outf.write("SRA_accession\tStatus\tOrganism\n")
        conn.close()


    def check_or_get_minikraken2(self, logger):
        """
        """
        if self.krakendir is None:
            self.krakendir = os.path.join(self.dbdir, "minikraken2_v2_8GB_201904_UPDATE", "")
        if not os.path.exists(self.krakendir):
            cmds = []
            if not os.path.exists(self.dbdir + "mini.tar.gz"):
                cmds.append(
                    str(
                        "wget ftp://ftp.ccb.jhu.edu/pub/data/kraken2_dbs/" +
                        "minikraken2_v2_8GB_201904_UPDATE.tgz -O {}mini.tar.gz"
                    ).format(self.dbdir))
            cmds.append("tar xzf {0}mini.tar.gz -C {0}".format(self.dbdir))
            logger.info("Downloading and preparing minikraken2 DB")
            for cmd in cmds:
                logger.debug(cmd)
                subprocess.run(
                    cmd,
                    shell=sys.platform != "win32",
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=True)


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
        conn = sqlite3.connect(self.db)
        c = conn.cursor()
        try:
            for acc, status, genus, species  in c.execute('SELECT * FROM SRAs'):
                # with open(self.SRAs_manifest, "r") as inf:
                # for i, line in enumerate(inf):
                # acc, status, org = line.strip().split("\t")
                self.SRAs[acc] = {
                    "status": status,
                    "genus": genus,
                    "species": species
                }
        except Exception as e:
            raise e
        conn.close()

    def update_manifest(self, newacc, newstatus, organism, logger):
        if " " in organism:
            genus, species = organism.split(" ")
        else:
            genus, species = organism, ""
        conn = sqlite3.connect(self.db)
        c = conn.cursor()

        #tmp = self.SRAs_manifest + ".bak"
        #shutil.move(self.SRAs_manifest, tmp)
        vals = (newacc, newstatus, genus, species, )
        c.execute('INSERT OR REPLACE INTO SRAs VALUES (?, ?, ?, ?)', vals)
        # with open(tmp, "r") as inf, open(self.SRAs_manifest, "w") as outf:
        #     for i, line in enumerate(inf):
        #         acc, status, lineorg = line.strip().split("\t")
        #         if acc == newacc:
        #             pass
        #         else:
        #             outf.write("{}\t{}\t{}\n".format(acc, status, lineorg))
            # if we still haven't written out our new SRA (ie, if we are adding
        # a new one, not updating)
            # outf.write("{}\t{}\t{}\n".format(newacc, newstatus, organism))
        # try:
        #     os.remove(tmp)
        # except FileNotFoundError:
        #     logger.warning("Missing backup manifest; multiple processes " +
                           # " could be trying to update it.")
        conn.commit()
        conn.close()
        self.read_SRA_manifest()

    def fetch_sraFind_data(self, logger):
        if self.sraFind_data is None:
            self.sraFind_data = os.path.join(
                self.dbdir, "sraFind.tab")
        sraFind_results = str(
            "https://raw.githubusercontent.com/nickp60/sraFind/master/" +
            "sraFind.tab"
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

    def get_SRA_data(self, SRA, org, logger, timeout, process_partial,
                     retry_partial, tool="fasterq-dump"):
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
        assert tool in ["fastq-dump", "fasterq-dump"], \
            "unrecognized download tool"
        suboutput_dir_raw = os.path.join(self.dbdir, SRA, "")
        logfile = os.path.join(suboutput_dir_raw, "download.log")
        if retry_partial and SRA in self.SRAs.keys():
            if self.SRAs[SRA]["status"] == "PARTIAL DOWNLOAD":
                if os.path.exists(suboutput_dir_raw):
                    for f in glob.glob(os.path.join(suboutput_dir_raw, "*.fastq")):
                        os.remove(f)
        rawreadsf, rawreadsr, download_error_message = \
            self.check_fastq_dir(this_data=suboutput_dir_raw,
                                 mate_as_single=False, logger=logger)
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
                    for f in glob.glob(os.path.join(suboutput_dir_raw, "*.fastq")):
                        os.remove(f)
            elif self.SRAs[SRA]['status'] == "LIBRARY TYPE ERROR":
                # dont try to reprocess
                return (None, None, "Library type Error")
        elif os.path.exists(suboutput_dir_raw):
            for f in glob.glob(os.path.join(suboutput_dir_raw, "*.fastq")):
                os.remove(f)
        else:
            pass
        os.makedirs(suboutput_dir_raw, exist_ok=True)
        # defaults to 6 threads or whatever is convenient;
        # we suspect I/O limits using more in most cases,
        # so we don't give the user the option to increase this
        # https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump
        verb = "" if tool == "fastq-dump" else "-vvv "
        cmd = str(
            "{tool} -O " +
            "{suboutput_dir_raw} --split-files {verb}{SRA} > {logfile}"
        ).format(**locals())
        logger.info("Downloading %s", SRA)
        logger.debug(cmd)
        try:
            subprocess.run(cmd,
                           shell=sys.platform != "win32",
                           timeout=timeout,
                           check=True)
        except subprocess.CalledProcessError:
            self.update_manifest(
                newacc=SRA,
                newstatus="DOWNLOAD ERROR",
                organism=org,
                logger=logger)
            logger.critical("Error running fasterq-dump; see log file at %s",
                            logfile)
            raise fasterqdumpError
        except subprocess.TimeoutExpired:
            if not process_partial:
                self.update_manifest(
                    newacc=SRA,
                    newstatus="DOWNLOAD ERROR",
                    organism=org,
                    logger=logger)
                logger.critical("%s timed out downloading %s",
                                tool, SRA)
                # delete any partial files;  if we try this right away,
                # fastq-dump doesn;t dump out fast enough, and we have nothing to
                # delete, but files apear later.  So we sleep for a few
                logger.debug("removing incomplete files from fast(er)q-dump")
                time.sleep(10)
                for f in glob.glob(suboutput_dir_raw + "*.fastq"):
                    os.remove(f)
                raise fasterqdumpError
            else:
                logger.warning("Downloaded %s partially; continuing", SRA)
                time.sleep(10)
                self.update_manifest(
                    newacc=SRA,
                    newstatus="PARTIAL DOWNLOAD",
                    organism=org,
                    logger=logger)
        rawreadsf, rawreadsr, download_error_message = \
            self.check_fastq_dir(this_data=suboutput_dir_raw,
                                 mate_as_single=False, logger=logger)
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

    def check_fastq_dir(self, this_data, mate_as_single, logger):
        # Double check the download worked.  If its a single lib,
        # it wont have a _1 prefix, so we check if that exists
        # and if so, adjust expeactations
        if not os.path.exists(this_data):
            # this never happens unless a run is aborted;
            # regardless, we want to make sure we attempt to re-download
            return(None, None, "No directory created for SRA during download")
        fastqs = glob.glob(os.path.join(this_data, "", "*.f*"))
        logger.debug("fastqs detected: %s", " ".join(fastqs))
        if len(fastqs) == 0:
            return(None, None, "No fastq files downloaded")
        if fastqs[0].endswith(".fastq"):
            pre = ".fastq"
        elif fastqs[0].endswith(".gz"):
            pre = ".fastq.gz"
        else:
            return(
                None, None,
                "Unfamiliar prefix for %s: must be either fastq or .gz" %
                fastqs[0])
        rawf, rawr = [], []
        rawreadsf, rawreadsr = None, None
        download_error_message = ""
        for fastq in fastqs:
            if fastq.endswith("_1" + pre):
                # not appending, cause we want to bump out any single libs t
                # that may have been read in first
                if len(rawf) != 0:
                    logger.warning("ignoring extra library %s", " ".join(rawf))
                rawf = [fastq]
            elif fastq.endswith("_2" + pre):
                rawr.append(fastq)
            elif fastq.endswith(pre) and not fastq.endswith("_3" + pre):
                if len(rawf) == 0:
                    # This is how we treat single libraries
                    rawf.append(fastq)
                else:
                    logger.warning("ignoring extra library %s", fastq)
            else:
                logger.error("Can only process paired or single-end reads")
                if mate_as_single:
                    logger.warning(
                        "treating first read from mate-pair library as single")
                else:
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
            if not rawreadsf.endswith("_1" + pre) and rawreadsr is not None:
                download_error_message = "cannot process a single library " + \
                    " (non-forward) file and a reverse file"
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
        ngenomes = len(glob.glob(os.path.join(self.refdir, "*.fna")))
        # check how many we need - we need to have at least the number
        #  in args.nstrains. If some exist, we re-rerun with the difference
        if ngenomes == 0 or args.nstrains > ngenomes:
            logger.info('Downloading genomes')
            return(self.our_get_n_genomes(
                org=args.organism_name,
                nstrains=args.nstrains - ngenomes,
                thisseed=args.seed,
                logger=logger)
            )
        # if len(glob.glob(os.path.join(self.refdir, "*.fna.gz"))) != 0:
        #     logger.warning('Genome downloading may have been interupted; ' +
        #                    'downloading fresh')
        #     assert self.refdir != self.dbdir, "critical: refdir and dbdir should never be the same!"
        #     shutil.rmtree(self.refdir)
        #     os.makedirs(self.refdir)
        #     return(self.our_get_n_genomes(
        #         org=args.organism_name,
        #         nstrains=args.nstrains,
        #         thisseed=args.seed,
        #         logger=logger)
        #     )
        return 0

    def our_get_n_genomes(self, org, nstrains,  thisseed, logger):
        # taken from the main function of get_n_genomes
        if self.prokaryotes is None:
            self.prokaryotes = os.path.join(self.dbdir, "prokaryotes.txt")
        if not os.path.exists(self.prokaryotes):
            gng.fetch_prokaryotes(dest=self.prokaryotes)
        org_lines = gng.get_lines_of_interest_from_proks(
            path=self.prokaryotes, org=org)
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
        # this can happen if tabs end up in metadata (see
        # AP019724.1 B. unifomis, and others
        # I updated plentyofbugs to try to catch it, but still might happen
        if len(cmds) == 0:
            return 1
        fetched = 0
        try:
            for i, cmd in enumerate(cmds):
                # ignore if already there; the last par of the  command is the zipped genome:
                thisgenome = cmd.split(" ")[-1].replace(".gz", "")
                if os.path.exists(os.path.join(self.refdir,  thisgenome)):
                    logger.debug(thisgenome + " already present, skipping")
                    continue
                else:
                    fetched = fetched + 1
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
        if fetched == 0:
            return 0
        for gz in glob.glob(os.path.join(self.refdir, "*.gz")):
            unzip_cmd = "gunzip " + gz
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
