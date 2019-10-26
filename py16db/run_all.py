#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import subprocess
import shutil
import gzip
import logging
import glob
import multiprocessing

from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from py16db.run_sickle import run_sickle

from . import __version__
from py16db.FocusDBData import FocusDBData, fasterqdumpError
from py16db.shared_methods import filter_SRA


class bestreferenceError(Exception):
    pass


class coverageError(Exception):
    pass


class kraken2Error(Exception):
    pass


class referenceNotGoodEnoughError(Exception):
    pass


class downsamplingError(Exception):
    pass


class riboSeedError(Exception):
    pass


class riboSeedUnsuccessfulError(Exception):
    """ its not magic, this "error" is for when riboSeed
    finishes, but cant improve on assembly
    """
    pass


class extracting16sError(Exception):
    pass


class barrnapError(Exception):
    pass


class libraryError(Exception):
    pass


def setup_logging(args):  # pragma: nocover
    if (args.verbosity * 10) not in range(10, 60, 10):
        raise ValueError('Invalid log level: %s' % args.verbosity)
    logging.basicConfig(
        level=logging.DEBUG,
        filemode='w',
        filename=os.path.join(args.output_dir, "focusDB.log"),
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    logger = logging.getLogger(__name__)
    console_err = logging.StreamHandler(sys.stderr)
    console_err.setLevel(level=(args.verbosity * 10))
    console_err_format = logging.Formatter(
        str("%(asctime)s \u001b[3%(levelname)s\033[1;0m  %(message)s"),
        "%H:%M:%S")
    console_err.setFormatter(console_err_format)
    logging.addLevelName(logging.DEBUG,    "4m --")
    logging.addLevelName(logging.INFO,     "2m ==")
    logging.addLevelName(logging.WARNING,  "3m !!")
    logging.addLevelName(logging.ERROR,    "1m xx")
    logging.addLevelName(logging.CRITICAL, "1m XX")
    logger.addHandler(console_err)
    return logger


def get_args():  # pragma: nocover
    parser = argparse.ArgumentParser(
        description="For a given genus or species, " +
        "focusDB orchestrates downloading whole-genome sequencing SRA, " +
        "running quality conrol, taxonimic assignment, reassembling with " +
        "riboSeed and extraction of 16s sequences",
        add_help=False)  # to allow for custom help
    mainargs = parser.add_argument_group('Main options')
    parargs = parser.add_argument_group('Parameters')
    jobargs = parser.add_argument_group('Job Handling')
    configargs = parser.add_argument_group('Configuration')
    expargs = parser.add_argument_group('Expert')
    #  these args are solely used internally
    hiddenargs = parser.add_argument_group('Hidden')
    mainargs.add_argument(
        "-o", "--output_dir",
        help="path to output", required=True)
    mainargs.add_argument(
        "-", "--organism_name",
        help="genus or genus species in quotes",
        required=True)

    hiddenargs.add_argument(
        "-g", "--genus",
        help=argparse.SUPPRESS,
        required=False)
    hiddenargs.add_argument(
        "-s", "--species",
        help=argparse.SUPPRESS,
        required=False)
    expargs.add_argument(
        "--SRA_list",
        help="path to file containing list of sras " +
        "for assembly [one column]",
        required=False)
    mainargs.add_argument(
        "--SRAs", default=None, nargs="+",
        help="negates -n_SRAs; " +
        "run pipeline on this (these) SRA(s) only",
        required=False)
    mainargs.add_argument(
        "-h", "--help",
        action="help", default=argparse.SUPPRESS,
        help="Displays this help message")

    parargs.add_argument(
        "-S", "--n_SRAs", help="max number of SRAs to be run",
        type=int, required=False)
    parargs.add_argument(
        "-R", "--n_references",
        help="max number of reference strains to consider. " +
        "default (0) is download all",
        type=int, required=False, default=0)

    configargs.add_argument(
        "--version", action='version',
        version='%(prog)s {version}'.format(
            version=__version__))
    parargs.add_argument(
        "-l", "--approx_length",
        help="Integer for approximate genome length",
        required=False, type=int)
    parargs.add_argument(
        "--just_seed", action="store_true",
        help="run just the subassemblies; " +
        "shorter execution time",
        required=False)
    configargs.add_argument(
        "--sraFind_path", dest="sra_path",
        help="path to sraFind file; " +
        "default is in ~/.focusDB/",
        default=None,
        required=False)
    configargs.add_argument(
        "--focusDB_data", dest="focusDB_data",
        help="path to data storage area; default ~/.focusDB/",
        default=None)
    configargs.add_argument(
        "--genomes_dir",
        help="path to where reference genomes are/will be " +
        "stored . Default location " +
        "is ~/.focusDB/references/genus_species/")
    #  Note  this arg doesn't get called, but is inheirited by get_n_genomes
    configargs.add_argument(
        "--prokaryotes", action="store",
        help="path to prokaryotes.txt; default is " +
        "in ~/.focusDB/",
        default=None,
        required=False)
    configargs.add_argument(
        "--kraken2_dir", action="store",
        help="path to kraken dir; default is " +
        "in ~/.focusDB/.  Will be created if doesn't exist",
        default=None,
        required=False)
    configargs.add_argument(
        "--kraken_mem_mapping", action="store_true",
        help="use this flag to load kraken2 db via disk " +
        "instead of RAM for taxonomic assignment. " +
        "automatically enabled if --memory < 20GB",
        required=False)
    expargs.add_argument(
        "--get_all",
        help="if a biosample is associated with " +
        "multiple libraries, default behaviour is to " +
        "download the first only.  Use --get_all to " +
        "analyse each library",
        action="store_true", required=False)
    parargs.add_argument(
        "--maxdist",
        help="maximum mash distance allowed for reference " +
        "genome; defaults to .05 (see Mash paper), which " +
        "roughly corresponds to species level similarity. " +
        "If desired, this can be relaxed",
        default=.05,
        type=float)
    jobargs.add_argument(
        "--njobs",
        help="how many jobs to run concurrently " +
        "via multiprocessing. --cores and --memory is per job",
        default=1, type=int)
    jobargs.add_argument(
        "--cores",
        help="PER JOB: how many cores you wish to use",
        default=1,
        required=False, type=int)
    jobargs.add_argument(
        "--memory",
        help="PER JOB: amount of RAM to be used. riboSeed " +
        "needs 10GB ram to run optimally; less, riboSeed " +
        "runs in --serialize mode to prevent memory errors" +
        "during subassemblies",
        default=4,
        required=False, type=int)
    jobargs.add_argument(
        "--timeout",
        help="Download SRAs can stall out periodically; " +
        "it typically takes 5-15 minutes for an average SRA" +
        "default 1800s (30 mins)",
        default=1800,
        required=False, type=int)
    expargs.add_argument(
        "--process_partial",
        help="If fastq-dump (NOT fasterq-dump) times out, " +
        "process what has been downloaded so far. This is " +
        "useful when timeout exceeds cause an SRA is very " +
        "large.  Becasue we will downsample later, enabling " +
        "this option allows processing of the partial " +
        "file(s). Default is to delete the partial files " +
        "and retry next time. Consider increasing --mincov  " +
        "to ensure that you only process partial files of " +
        "sensible size. EXPERTS ONLY",
        required=False, action="store_true")
    expargs.add_argument(
        "--retry_partial",
        help="If a partial download is encountered during " +
        "this run, delete and attempt to re-download",
        required=False, action="store_true")
    jobargs.add_argument(
        "--threads",
        action="store",
        default=1, type=int,
        choices=[1, 2, 4],
        help="if your cores are hyperthreaded, set number" +
        " threads to the number of threads per processer." +
        "If unsure, see 'cat /proc/cpuinfo' under 'cpu " +
        "cores', or 'lscpu' under 'Thread(s) per core'." +
        ": %(default)s")
    parargs.add_argument(
        "--maxcov",
        help="integer for maximum desired read depth" +
        "after downsampling",
        default=50,
        required=False, type=int)
    parargs.add_argument(
        "--mincov",
        help="integer for minimum  read depth",
        default=15,
        required=False, type=int)
    jobargs.add_argument(
        "--fastqtool",
        help="either fastq-dump or fasterq-dump",
        default="fasterq-dump",
        choices=["fastq-dump", "fasterq-dump"],
        required=False)
    expargs.add_argument(
        "--custom_reads",
        help="input of custom reads", nargs='+',
        required=False, type=str)
    expargs.add_argument(
        "--custom_name",
        help="if using --custom_reads, store as this name",
        required=False, type=str)
    expargs.add_argument(
        "--redo_assembly", action="store_true",
        help="redo the assembly step, ignoring status file")
    jobargs.add_argument(
        "--subassembler",
        help="which program should riboseed " +
        "use for sub assemblies",
        choices=["spades", "skesa"],
        required=False, default="spades")
    # this is needed for plentyofbugs, should not be user set
    hiddenargs.add_argument(
        "--nstrains", help=argparse.SUPPRESS,
        type=int, required=False)
    expargs.add_argument(
        "--seed",
        help="random seed for subsampling references and SRAs",
        type=int, default=12345)
    jobargs.add_argument(
        "-v", "--verbosity", dest='verbosity',
        action="store",
        default=2, type=int, choices=[1, 2, 3, 4, 5],
        help="Logger writes debug to file in output dir; " +
        "this sets verbosity level sent to stderr. " +
        " 1 = debug(), 2 = info(), 3 = warning(), " +
        "4 = error() and 5 = critical(); " +
        "default: %(default)s")
    args = parser.parse_args()
    if " " in args.organism_name:
        args.genus, args.species = args.organism_name.split(" ")
    else:
        args.genus, args.species = args.organism_name, ""
    if args.custom_reads is not None:
        if args.custom_name is None:
            print("--custom_name is required using custom reads")
            print("this name is used to store the reads in the focusDB data ")
            sys.exit(1)
        else:
            args.custom_name = args.custom_name.replace(" ", "_")
    # plentyofbugs uses args.nstrains, but we call
    # it args.n_references for clarity
    args.nstrains = args.n_references
    if args.SRAs is None:
        if args.custom_reads is None:
            if args.n_SRAs is None:
                print("if not running with --SRAs, " +
                      "then --n_SRAs must be provided!")
                sys.exit(1)
    # catch process partial and fasterq-dump
    if args.process_partial and args.fastqtool == "fasterq-dump":
        print("--process_partial can only be used with --fastqtool fastq-dump")
        sys.exit(1)
    return(args)


def check_programs(logger):
    """exits if the following programs are not installed"""

    required_programs = [
        "ribo", "barrnap", "fasterq-dump", "mash",
        "skesa", "plentyofbugs", "iqtree", "seqtk",
        "kraken2"]
    for program in required_programs:
        if shutil.which(program) is None:
            logger.critical('%s is not installed: exiting.', program)
            sys.exit(1)


def parse_status_file(path):
    # because downloading and assembling can fail for many reasons,
    # we write out status to a file.  this allows for easier restarting of
    # incomplete runs
    if not os.path.exists(path):
        return []
    statuses = []
    with open(path, "r") as statusfile:
        for line in statusfile:
            statuses.append(line.strip())
    return(statuses)


def update_status_file(path, to_remove=[], message=None):
    assert isinstance(to_remove, list)
    statuses = parse_status_file(path)
    # dont try to remove empty files
    if statuses != []:
        os.remove(path)
    # just for cleaning up status file
    if message is not None:
        statuses.append(message)
    # write out non-duplicated statuses
    with open(path, "w") as statusfile:
        for status in set(statuses):
            if status not in to_remove:
                statusfile.write(status + "\n")


def sralist(list):
    """ takes a file list of  of SRAs, return list
    for if you wish to use SRAs that are very recent and ahead of sraFind
    """
    sras = []
    with open(list, "r") as infile:
        for sra in infile:
            sras.append(sra.strip())
    return sras


def pob(genomes_dir, readsf, output_dir, maxdist, logger):
    """use plentyofbugs to identify best reference
    Uses plentyofbugs, a package that useqs mash to
    find the best reference genome for draft genome
    """
    pobcmd = str("plentyofbugs -g {genomes_dir} -f {readsf} -o {output_dir} " +
                 "--downsampling_ammount 1000000").format(**locals())
    logger.info('Finding best reference genome: %s', pobcmd)

    for command in [pobcmd]:
        logger.debug(command)
        try:
            subprocess.run(command,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
            best_ref = os.path.join(output_dir, "best_reference")
        except Exception as e:
            logger.error(e)
            raise bestreferenceError(
                "Error running the following command: %s" % command)

    with open(best_ref, "r") as infile:
        for line in infile:
            sraacc = line.strip().split('\t')
            sim = float(sraacc[1])
            ref = sraacc[0]
            logger.debug("Reference genome mash distance: %s", sim)
            if sim > maxdist:
                raise referenceNotGoodEnoughError(
                    "Reference similarity %s does not meet %s threshold" %
                    (sim, maxdist))
            length_path = os.path.join(output_dir, "genome_length")
            cmd = "wc -c {ref} > {length_path}".format(**locals())
            subprocess.run(cmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
            return(ref, sim)


def check_rDNA_copy_number(ref, output, logger):
    """ensure reference has multiple rDNAs
    Using barrnap to check that there are multiple rDNA copies
    in the reference genome
    """
    os.makedirs(os.path.join(output, "barrnap_reference"), exist_ok=True)
    barroutput = os.path.join(output, "barrnap_reference",
                              os.path.basename(ref) + ".gff")
    cmd = "barrnap {ref} > {barroutput}".format(**locals())
    subprocess.run(cmd,
                   shell=sys.platform != "win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    rrn_num = 0
    with open(barroutput, "r") as rrn:
        for rawline in rrn:
            line = rawline.strip().split('\t')
            if line[0].startswith("##"):
                continue
            if line[8].startswith("Name=16S"):
                rrn_num += 1
    return rrn_num


def get_and_check_ave_read_len_from_fastq(fastq1, minlen, maxlen, logger=None):
    """return average read length in fastq1 file from first N reads
    from LP: taken from github.com/nickp60/riboSeed/riboSeed/classes.py;
    """
    nreads = None
    # nreads = 30
    tot = 0
    if os.path.splitext(fastq1)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    logger.debug("Obtaining average read length")
    with open_fun(fastq1, "rt") as file_handle:
        data = SeqIO.parse(file_handle, "fastq")
        for i, read in enumerate(data):
            tot += len(read)
            if nreads is not None:
                if i  >= nreads:
                    break

    ave_read_len = float(tot / i)
    if ave_read_len < minlen:
        logger.error("Average read length is too short: %s; skipping...",
                     ave_read_len)
        return (1, ave_read_len)
    if ave_read_len > maxlen:
        logger.critical("Average read length is too long: %s; skipping...",
                        ave_read_len)
        return (2, ave_read_len)
    logger.debug("Average read length: %s", ave_read_len)
    return (0, ave_read_len)


def get_coverage(read_length, approx_length, fastq1, fastq2, logger):
    """Obtains the coverage for a read set given the estimated genome size"""
    if os.path.splitext(fastq1)[-1] in ['.gz', '.gzip']:
        open_fun = gzip.open
    else:
        open_fun = open
    logger.debug("Counting reads")

    with open_fun(fastq1, "rt") as data:
        for count, line in enumerate(data):
            pass

    if fastq2 is not None:
        read_length = read_length * 2

    coverage = float((count * read_length) / (approx_length * 4))
    logger.info('Read coverage: %sx', round(coverage, 1))
    return(coverage)


def downsample(read_length, approx_length, fastq1, fastq2,
               mincoverage, maxcoverage, destination, logger, run):
    """downsample for optimal assembly
    Given the coverage from coverage(), downsamples the reads if over
    the max coverage set by args.maxcov. Default 50.
    """
    suboutput_dir_downsampled = destination
    downpath1 = os.path.join(suboutput_dir_downsampled,
                             "downsampledreadsf.fastq")
    downpath2 = None
    if fastq2 is not None:
        downpath2 = os.path.join(suboutput_dir_downsampled,
                                 "downsampledreadsr.fastq")
    return_originals = True
    if not run:
        # if any downsmapled reads are here, return those;
        # othewise, we assume they did not need to be downsmapled
        for f in [downpath1, downpath2]:
            if f is not None:
                if os.path.exists(f):
                    return_originals = False
        if return_originals:
            return(fastq1, fastq2)
        else:
            return(downpath1, downpath2)
    coverage = get_coverage(read_length, approx_length,
                            fastq1, fastq2, logger=logger)
    if coverage < mincoverage:
        raise coverageError("%sx coverage fails to meet minimum (%s)" %
                            (coverage, mincoverage))
    # seqtk either works with a number of reads, or a fractional value
    # for how many reads to retain.  Here we calculate the later based
    # on what we have currently
    covfrac = round(float(maxcoverage / coverage), 3)
    stk_cmd_p = "seqtk sample -s100"
    dcmd = "{stk_cmd_p} {fastq1} {covfrac} > {downpath1}".format(**locals())
    dcmd2 = "{stk_cmd_p} {fastq2} {covfrac} > {downpath2}".format(**locals())
    # at least downsample the forward/single reads, but add the
    # other command if using paired reads
    commands = [dcmd]
    if (coverage > maxcoverage):
        if run:
            os.makedirs(suboutput_dir_downsampled)
            logger.info('Downsampling to %s X coverage', maxcoverage)
            if fastq2 is not None:
                commands.append(dcmd2)
            for command in commands:
                try:
                    logger.debug(command)
                    subprocess.run(command,
                                   shell=sys.platform != "win32",
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE,
                                   check=True)
                except Exception as e:
                    logger.error(e)
                    raise downsamplingError(
                        "Error running following command ", command)
        return(downpath1, downpath2)
    else:
        logger.info(
            'Skipping downsampling as max coverage is < %s', maxcoverage)
        return(fastq1, fastq2)


def make_riboseed_cmd(sra, readsf, readsr, cores, subassembler, threads,
                      output, memory, just_seed, logger):
    """Runs riboSeed to reassemble reads """
    if memory < 10:
        serialize = "--serialize "
    else:
        serialize = ""
    if just_seed:
        seed_str = "--just_seed "
    else:
        seed_str = ""
    cmd = str("ribo run -r {sra} -F {readsf} -R {readsr} --cores {cores} " +
              "--threads {threads} -v 1 -o {output} {serialize}" +
              "--subassembler {subassembler} {seed_str}--skip_control " +
              "--stages none " +
              "--memory {memory}").format(**locals())

    if readsr is None:
        cmd = str("ribo run -r {sra} -S1 {readsf} --cores {cores} " +
                  "--threads {threads} -v 1 -o {output} {serialize}" +
                  "--subassembler {subassembler} {seed_str}--skip_control " +
                  "--stages none " +
                  "--memory {memory}").format(**locals())
    return(cmd)


def run_kraken2(args, contigs, dest_prefix, db, logger):
    out = dest_prefix + ".output"
    report = dest_prefix + ".report"
    if args.memory < 20 or args.kraken_mem_mapping:
        memstring = "--memory-mapping "
    else:
        memstring = ""
    cmd = str(
        "kraken2 {memstring}--db {db} --threads {args.cores} " +
        "--use-names --output {out} " +
        "--report {report} {contigs}").format(**locals())
    if not os.path.exists(report):
        logger.debug(cmd)
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    return report


def parse_kraken_report(kraken2_report):
    """column data:
    https://ccb.jhu.edu/software/kraken2/index.shtml?t=manual
    #standard-kraken-output-format
    """
    empty = [0, "-", ""]
    tax = {'R': empty,
           'D': empty,
           'P': empty,
           'C': empty,
           'O': empty,
           'F': empty,
           'G': empty,
           'S': empty}
    with open(kraken2_report, "r") as inf:
        for line in inf:
            # print(line)
            sline = line.split("\t")
            if len(sline) != 6:
                raise ValueError(
                    "Malformed kraken2 report; should be 6 columns: %s" % line)
            perc, n_in, n_at, lev, taxid, label = sline
            if lev in tax.keys():
                # only report the top hit
                this = float(perc.strip()), taxid, label.strip()
                if tax[lev][0] < this[0]:
                    tax[lev] = float(perc.strip()), taxid, label.strip()
    return tax


def process_strain(rawreadsf, rawreadsr, read_length, genomes_dir,
                   this_output, args, logger, status_file, kdir):
    """return a tuple of the riboSeed cmd and the path to contigs,
    and the taxonomy according to kraken
    """
    pob_dir = os.path.join(this_output, "plentyofbugs", "")
    krak_dir = os.path.join(this_output, "kraken2", "")
    ribo_dir = os.path.join(this_output, "riboSeed", "")
    sickle_out = os.path.join(this_output, "sickle",  "")
    best_reference = os.path.join(pob_dir, "best_reference")

    # Note thhat the status file is checked before each step.
    # If a failure occured, all future steps are rerrun
    # for instance, if trimming has been done, but downsample hasn't,
    # downsampling and assembly will be run. This is to protect against
    # files sticking around when they shouldn't
    if not os.path.exists(best_reference):
        if os.path.exists(pob_dir):
            shutil.rmtree(pob_dir)
        pob(genomes_dir=genomes_dir, readsf=rawreadsf,
            output_dir=pob_dir, maxdist=args.maxdist, logger=logger)

    with open(best_reference, "r") as infile:
        for line in infile:
            best_ref_fasta = line.split('\t')[0]
            best_ref_dist = float(line.split('\t')[1])
        if args.maxdist < best_ref_dist:
            raise referenceNotGoodEnoughError(
                "Reference similarity %s does not meet %s threshold" %
                (best_ref_dist, args.maxdist))
    report_output = krak_dir + "kraken2.report"
    if "TAXONOMY" not in parse_status_file(status_file) or \
       not os.path.exists(report_output):
        if os.path.exists(krak_dir):
            shutil.rmtree(krak_dir)
        os.makedirs(krak_dir, exist_ok=True)
        logger.info('Assigning taxonomy with kraken')
        pob_assembly = os.path.join(pob_dir, "assembly", "contigs.fasta")
        try:
            report_output = run_kraken2(
                args,
                contigs=pob_assembly,
                dest_prefix=krak_dir + "kraken2",
                db=kdir, logger=logger)
            update_status_file(status_file, message="TAXONOMY")
        except Exception as e:
            raise kraken2Error(e)
    if os.path.getsize(report_output) == 0:
        raise kraken2Error("Kraken output file  exists but is empty")
    tax_dict = parse_kraken_report(kraken2_report=report_output)
    logger.debug(tax_dict)

    if args.approx_length is None:
        genome_length = os.path.join(pob_dir, "genome_length")
        approx_length = None
        with open(genome_length, "r") as infile:
            for line in infile:
                approx_length = float(line.split()[0])
                logger.debug("Using genome length: %s", approx_length)
        if approx_length is None:
            raise ValueError("Error running plentyofbugs; " +
                             "database possibly outdated")
    else:
        approx_length = args.approx_length
    if "TRIMMED" not in parse_status_file(status_file):
        logger.info('Quality trimming reads')
        update_status_file(status_file,
                           to_remove=["DOWNSAMPLED", "RIBOSEED COMPLETE"])
        if os.path.exists(sickle_out):
            shutil.rmtree(sickle_out)
    trimmed_fastq1, trimmed_fastq2 = run_sickle(
        fastq1=rawreadsf,
        fastq2=rawreadsr,
        output_dir=sickle_out,
        run="TRIMMED" not in parse_status_file(status_file))
    update_status_file(status_file, message="TRIMMED")
    logger.debug('Quality trimmed f reads: %s', trimmed_fastq1)
    logger.debug('Quality trimmed r reads: %s', trimmed_fastq2)

    # if if os.path.getsize(trimmed_fastq1) == 0:
    #     raise libraryError("Error occured when trimming. This could be due " +
    #                        "to incorrect metadata about pairing. " +
    #                        "For more information, see Sickle results in " +
    #                        sickle_out)
    logger.debug('Downsampling reads')
    if "DOWNSAMPLED" not in parse_status_file(status_file):
        update_status_file(status_file, to_remove=["RIBOSEED COMPLETE"])
        if os.path.exists(os.path.join(this_output, "downsampled")):
            shutil.rmtree(os.path.join(this_output, "downsampled"))
    downsampledf, downsampledr = downsample(
        approx_length=approx_length,
        fastq1=trimmed_fastq1,
        fastq2=trimmed_fastq2,
        mincoverage=args.mincov,
        maxcoverage=args.maxcov,
        destination=os.path.join(this_output, "downsampled"),
        read_length=read_length,
        logger=logger,
        run="DOWNSAMPLED" not in parse_status_file(status_file))
    update_status_file(status_file, message="DOWNSAMPLED")
    logger.debug('Downsampled f reads: %s', downsampledf)
    logger.debug('Downsampled r reads: %s', downsampledr)
    riboseed_cmd = make_riboseed_cmd(sra=best_ref_fasta, readsf=downsampledf,
                                     readsr=downsampledr, cores=args.cores,
                                     memory=args.memory,
                                     subassembler=args.subassembler,
                                     threads=args.threads, output=ribo_dir,
                                     just_seed=args.just_seed,
                                     logger=logger)
    # do we want to redo the assembly?
    if args.redo_assembly:
        update_status_file(status_file, to_remove=["RIBOSEED COMPLETE"])
    # file that will contain riboseed contigs
    if args.just_seed:
        ribo_contigs = os.path.join(
            this_output, "riboSeed", "seed",
            "final_long_reads", "riboSeedContigs.fasta")
    else:
        ribo_contigs = os.path.join(
            this_output, "riboSeed", "seed",
            "final_de_fere_novo_assembly", "contigs.fasta")
    if "RIBOSEED COMPLETE" not in parse_status_file(status_file):
        if os.path.exists(ribo_dir):
            shutil.rmtree(ribo_dir)
        return(riboseed_cmd, ribo_contigs, tax_dict)
    else:
        logger.info("Skipping riboSeed")
        return (None, ribo_contigs, tax_dict)


def check_riboSeed_outcome(status_file, contigs):
    # check for the files to see if riboSeed completed
    if os.path.exists(contigs):
        update_status_file(status_file, message="RIBOSEED COMPLETE")
    else:
        this_output = os.path.dirname(
            os.path.dirname(os.path.dirname(contigs)))
        raise riboSeedUnsuccessfulError(str(
            "riboSeed completed but was not successful; " +
            "for details, see log file at %s") %
            os.path.join(this_output, "run_riboSeed.log"))


def run_barrnap(assembly,  results, logger):
    barrnap = "barrnap {assembly} > {results}".format(**locals())
    logger.debug('Identifying 16S sequences with barnap: %s', barrnap)
    try:
        subprocess.run(barrnap,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    except Exception as e:
        logger.error(e)
        raise barrnapError(
            "Error running the following command %s" % barrnap)


def extract_16s_from_assembly(assembly, gff, sra, output, output_summary,
                              args, singleline, tax_d, min_length, logger):
    tax_string = tax_d["S"][2]
    # if no label, gibe the next one
    if len(tax_string.replace(" ", "")) == 0:  # if genus only
        tax_string = tax_d["G"][2] + "sp."
        if len(tax_string.replace(" ", "").replace(".sp", "")) == 0:
            tax_string = tax_d["F"][2]
            if len(tax_string.replace(" ", "")) == 0:
                tax_string = tax_d["O"][2]
                if len(tax_string.replace(" ", "")) == 0:
                    tax_string = tax_d["C"][2]
                    if len(tax_string.replace(" ", "")) == 0:
                        tax_string = tax_d["P"][2]
    score_string = str(tax_d["D"][0])
    taxid_string = tax_d["D"][1]
    big_tax_string = tax_d["D"][2]
    for lev in ["P", "C", "O", "F", "G", "S"]:
        score_string = score_string + ";" + str(tax_d[lev][0])
        taxid_string = taxid_string + ";" + tax_d[lev][1]
        big_tax_string = big_tax_string + ";" + tax_d[lev][2]

    # struction of results16s:
    # [sra_#, chromosome, start, end, reverse complimented,
    #  big_tax_string, score_string, taxid_string, tax_string]
    results16s = {}
    nseqs = 0
    with open(gff, "r") as rrn, open(output, "a") as outf, \
         open(output_summary, "a") as outsum:
        rrn_num = 0
        for rawline in rrn:
            line = rawline.strip().split('\t')
            # need this: catches index errors
            if line[0].startswith("##"):
                pass
            elif line[8].startswith("Name=16S"):
                rrn_num = rrn_num + 1
                if line[6] == "-":
                    suffix = 'chromosome-RC@'
                else:
                    suffix = ''
                chrom = line[0]
                ori = line[6]
                start = int(line[3])
                end = int(line[4])
                if (end - start) < 1000:
                    logger.info("ignoring short rRNA:")
                    logger.info(line)
                    continue
                thisid = "{}_{}".format(sra, rrn_num)
                results16s[thisid] = [chrom, start, end, line[6]]
                with open(assembly, "r") as asmb:
                    for rec in SeqIO.parse(asmb, "fasta"):
                        if rec.id != chrom:
                            continue
                        seq = rec.seq[start + 1: end + 1]
                        if ori == "-":
                            seq = seq.reverse_complement()
                        thisidcoords = "{thisid}.{start}.{end}".format(
                            **locals())
                        # Need to disable linewrapping for use with SILVA, etc
                        if singleline:
                            seqstr = str(seq)
                            outf.write(
                                str(">{thisidcoords} {tax_string}\n" +
                                    "{seqstr}\n").format(**locals()))
                        else:
                            SeqIO.write(
                                SeqRecord(
                                    seq, id=thisidcoords,
                                    description=tax_string),
                                outf,  "fasta")
                        outsum.write(
                            str(
                                "{thisid}\t{assembly}\t" +
                                "{chrom}\t{start}\t{end}\t{big_tax_string}\t" +
                                "{score_string}\t{taxid_string}\t" +
                                "{tax_string}\n"
                            ).format(**locals()))
                        nseqs = nseqs + 1
    return nseqs


def write_pass_fail(args, stage, status, note):
    """
    format fail messages in tabular fomat:
    organism\tstage\tmessage
    """
    path = os.path.join(args.output_dir, "SUMMARY")
    org = args.organism_name
    with open(path, "a") as failfile:
        failfile.write("{}\t{}\t{}\t{}\n".format(org, status, stage, note))


def run_riboSeed_catch_errors(cmd, acc=None, args=None, status_file=None,
                              riboSeed_jobs=None):
    if cmd is None:
        for j in riboSeed_jobs:
            if j[0] == acc:
                j[4] = 0
        return 0
    try:
        subprocess.run(cmd,
                       shell=sys.platform != "win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    except subprocess.CalledProcessError:
        for j in riboSeed_jobs:
            if j[0] == acc:
                j[4] = 1

        write_pass_fail(args, status="FAIL",
                        stage=acc,
                        note="Unknown failure running riboSeed")
    for j in riboSeed_jobs:
        if j[0] == acc:
            j[4] = 0
    return 0


def write_this_config(args, this_config_file):
    args_to_write = ["maxdist", "subassembler", "maxcov", "just_seed"]
    argd = vars(args)
    with open(this_config_file, "w") as outf:
        for arg in args_to_write:
            outf.write("{}:{}\n".format(arg, argd[arg]))


def different_args(args, this_config_file, logger):
    """ Returns empty list if no args differ
    """
    diff_args = []
    args_to_write = ["maxdist", "subassembler", "maxcov", "just_seed"]
    old_config_dict = {}
    this_config_dict = vars(args)
    if not os.path.exists(this_config_file):
        raise ValueError(
            "No previous config file found %s; rerunning" % this_config_file)
    with open(this_config_file, "r") as f:
        for line in f:
            (key, val) = line.strip().split(":")
            old_config_dict[key] = val
    if len(old_config_dict) == 0:
        raise ValueError("Old config file empty; rerunning")
    for arg in args_to_write:
        # note that reading from the file makes all old args strings, so we
        # accomodate that
        try:
            if str(this_config_dict[arg]) != old_config_dict[arg]:
                logger.info(
                    "New parameter value for " +
                    "{} ({}) doesn't match old value ({})".format(
                        arg, this_config_dict[arg], old_config_dict[arg]))
                diff_args.append(arg)
        except KeyError:
            logger.info("New parameter value for --just_seed")
            # this can happen after updates -- play it safe and redo it
            diff_args.append(arg)
    return diff_args


def add_key_or_increment(d, k):
    if k in d.keys():
        d[k] = d[k] + 1
    else:
        d[k] = 1


def main():
    args = get_args()

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    if os.path.exists(os.path.join(args.output_dir, "SUMMARY")):
        os.remove(os.path.join(args.output_dir, "SUMMARY"))

    logger = setup_logging(args)
    logger.info("Processing %s", args.organism_name)
    logger.info("Usage:\n{0}\n".format(" ".join([x for x in sys.argv])))
    logger.debug("All settings used:")
    for k, v in sorted(vars(args).items()):
        logger.debug("{0}: {1}".format(k, v))
    check_programs(logger)
    # set up the data object
    # grooms path names or uses default location if unset
    fDB = FocusDBData(
        dbdir=args.focusDB_data,
        refdir=args.genomes_dir,
        sraFind_data=args.sra_path,
        krakendir=args.kraken2_dir,
        prokaryotes=args.prokaryotes)
    fDB.check_genomes_dir(org=args.organism_name)
    fDB.fetch_sraFind_data(logger=logger)
    fDB.check_or_get_minikraken2(logger=logger)

    # process data 1 of 4 ways: specific SRA(s), a file of SRA(s),
    #  specific read file (stored as a faux SRA), or the default
    #  to get a list of SRAs from sraFind for a given organism name
    if args.SRAs is not None:
        filtered_sras = args.SRAs
    elif args.SRA_list is not None:
        filtered_sras = sralist(list=args.SRA_list)
    elif args.custom_reads is not None:
        this_data_dir = os.path.join(fDB.dbdir, args.custom_name)
        if not os.path.exists(this_data_dir):
            os.makedirs(this_data_dir)
            for read in args.custom_reads:
                dest = os.path.join(this_data_dir, os.path.basename(read))
                shutil.copy(read, dest)
        else:
            pass
        filtered_sras = [args.custom_name]
    else:
        filtered_sras = filter_SRA(
            sraFind=fDB.sraFind_data,
            organism_name=args.organism_name,
            strains=args.n_SRAs,
            thisseed=args.seed,
            logger=logger,
            get_all=args.get_all)

    if filtered_sras == []:
        if args.custom_reads is None:
            logger.critical('No SRAs found on NCBI by sraFind')
            write_pass_fail(
                args, status="FAIL",
                stage="global",
                note="No SRAs available")
            sys.exit(1)

    pob_result = fDB.decide_skip_or_download_genomes(args, logger)
    if pob_result != 0:
        if pob_result == 1:
            message = "No available references"
        elif pob_result == 2:
            message = "Error downloading genome from NCBI"
        elif pob_result == 3:
            message = "Error unzipping genomes; delete directory and try again"
        else:
            pass
        logger.critical(message)
        write_pass_fail(args, status="ERROR", stage="global", note=message)
        sys.exit(1)
    genome_check_file = os.path.join(fDB.refdir, ".references_passed_checks")
    # ########  Check to see if we have requested a different number of strains
    this_config_file = os.path.join(args.output_dir, "config")
    try:
        updated_args = different_args(args, this_config_file, logger)
    except ValueError as e:
        logger.warning(e)
        # if we have any issues finding or reading the config, rerun it all
        updated_args = ["maxdist", "subassembler", "maxcov", "just_seed"]
    # if "n_references" in updated_args:
    #     if os.path.exists(genome_check_file):
    #         os.remove(genome_check_file)
    write_this_config(args, this_config_file)

    # #########
    if not os.path.exists(genome_check_file):
        logger.info("checking reference genomes for rDNA counts")
        for pot_reference in glob.glob(os.path.join(fDB.refdir, "*.fna")):
            rDNAs = check_rDNA_copy_number(ref=pot_reference,
                                           output=fDB.refdir,
                                           logger=logger)
            if rDNAs < 2:
                logger.warning(
                    "reference %s does not have multiple rDNAs; excluding",
                    pot_reference)
                os.remove(pot_reference)
        with open(genome_check_file, "w") as statusfile:
            statusfile.write("References have been checked\n")
    else:
        logger.debug("Already checked reference genomes in %s", fDB.refdir)
    if len(glob.glob(os.path.join(fDB.refdir, "*.fna"))) == 0:
        logger.critical("No usable reference genome found!")
        write_pass_fail(args, status="FAIL",
                        stage="global",
                        note="No references had more than 1 rDNA")
        sys.exit(0)

    riboSeed_jobs = []  # [accession, cmd, contigs, status_file, return_code]
    nsras = len(filtered_sras)
    n_errors = {}
    for i, accession in enumerate(filtered_sras):
        # trying to troublshoot a potential race condition
        # deleting all references.
        assert len(glob.glob(os.path.join(fDB.refdir, "*.fna"))) != 0, \
            "as of SRA %s (%i of %i), genomes dir empty" % (
                accession, i + 1, nsras)
        this_output = os.path.join(args.output_dir, accession)
        this_results = os.path.join(this_output, "results")
        os.makedirs(this_output, exist_ok=True)
        status_file = os.path.join(this_output, "status")
        logger.info("Organism: %s; Accession: %s (%s of %s)",
                    args.organism_name, accession, i + 1, nsras)
        message = ""
        # ############### check updated args, update status file if needed
        if "maxdist" in updated_args:
            # plentyofbugs will rerun if this fileis missing
            update_status_file(status_file, to_remove=["RIBOSEED COMPLETE"])
            this_pob_results = os.path.join(
                this_results, "plentyofbugs", "best_reference")
            if os.path.exists(this_pob_results):
                os.remove(this_pob_results)
        if "maxcov" in updated_args:
            update_status_file(status_file,
                               to_remove=["DOWNSAMPLED", "RIBOSEED COMPLETE"])
        if "subassembler" in updated_args:
            update_status_file(status_file, to_remove=["RIBOSEED COMPLETE"])
        if "just_seed" in updated_args:
            update_status_file(status_file, to_remove=["RIBOSEED COMPLETE"])

        ################
        if "RIBOSEED COMPLETE" in parse_status_file(status_file) and \
           not args.redo_assembly:
            logger.info("using existing results")
            if args.just_seed:
                ribo_contigs = os.path.join(
                    this_results, "riboSeed", "seed",
                    "final_long_reads", "riboSeedContigs.fasta")
            else:
                ribo_contigs = os.path.join(
                    this_results, "riboSeed", "seed",
                    "final_de_fere_novo_assembly", "contigs.fasta")
            kraken2_report_output = os.path.join(
                this_results, "kraken2", "kraken2.report")
            # double check files exist before we skip this one
            if all([os.path.exists(p) and os.path.getsize(p) != 0
                    for p in [ribo_contigs, kraken2_report_output]]):
                riboSeed_jobs.append(
                    [accession, None, ribo_contigs,  status_file,
                     parse_kraken_report(kraken2_report_output), 0])
                continue

        try:
            rawreadsf, rawreadsr, download_error_message = \
                fDB.get_SRA_data(
                    org=args.organism_name,
                    # genus=args.genus,    # TODO if needed
                    # species=args.species
                    SRA=accession,
                    logger=logger,
                    timeout=args.timeout,
                    process_partial=args.process_partial,
                    retry_partial=args.retry_partial,
                    tool=args.fastqtool)
        except fasterqdumpError:
            message = 'Error downloading %s' % accession
            write_pass_fail(
                args, status="ERROR", stage=accession, note=message)
            logger.error(message)
            add_key_or_increment(n_errors, "Downloading")
            continue
        if download_error_message != "":
            write_pass_fail(args, status="ERROR", stage=accession,
                            note=download_error_message)
            logger.error(
                "Error either downloading or parsing the file " +
                "name for this accession.")
            logger.error(download_error_message)
            continue
        read_len_status, read_length = get_and_check_ave_read_len_from_fastq(
            minlen=65,
            maxlen=301,
            fastq1=rawreadsf, logger=logger)
        if read_len_status != 0:
            if read_len_status == 1:
                message = "Reads were shorter than 65bp threshold"
            else:
                message = "Reads were longer than 300bp threshold"
            write_pass_fail(
                args, status="ERROR", stage=accession, note=message)
            logger.error(message)
            continue
        #  heres the meat of the main, catching errors for
        #  anything but the riboSeed step
        try:
            riboSeed_cmd, contigs_path, taxonomy_d = process_strain(
                rawreadsf, rawreadsr, read_length, fDB.refdir,
                this_results, args, logger, status_file, fDB.krakendir)
            riboSeed_jobs.append(
                [accession, riboSeed_cmd,
                 contigs_path,  status_file, taxonomy_d, None])
        except coverageError as e:
            write_pass_fail(args, status="FAIL",
                            stage=accession,
                            note="Insufficient coverage")
            logger.error(e)
            continue
        except bestreferenceError as e:
            write_pass_fail(args, status="ERROR",
                            stage=accession,
                            note="Unknown error selecting reference")
            logger.error(e)
            add_key_or_increment(n_errors, "plentyofbugs")
            continue
        except kraken2Error as e:
            if not args.kraken_mem_mapping:
                logger.error("Kraken2 error; try rerunning with " +
                             "--kraken_mem_mapping")
            write_pass_fail(args, status="ERROR",
                            stage=accession,
                            note="Unknown error runing kraken2")
            add_key_or_increment(n_errors, "Taxonomy")
            logger.error(e)
            continue
        except referenceNotGoodEnoughError as e:
            write_pass_fail(
                args, status="FAIL",
                stage=accession,
                note="No reference meets threshold for re-assembly")
            logger.error(e)
            continue
        except downsamplingError as e:
            write_pass_fail(args, status="ERROR",
                            stage=accession,
                            note="Unknown error downsampling")
            logger.error(e)
            add_key_or_increment(n_errors, "Downsampling")
            continue
        except Exception as e:
            logger.error(e)
            logger.error(
                "Unknown error occured; please raise issue on GitHub " +
                "attaching the log file found in %s ", this_results
                )
            add_key_or_increment(n_errors, "Unknown")
            write_pass_fail(args, status="FAIL",
                            stage=accession,
                            note="Unknown critial error")
            continue

    #######################################################################
    all_assemblies = []  # [contigs, tax{}]
    ribo_cmds = [x[1] for x in riboSeed_jobs if x[1] is not None]
    # split_cores = int(args.cores / (len(ribo_cmds) / 2))
    # if split_cores < 1:
    #     split_cores = 1
    if len(riboSeed_jobs) > 0:
        logger.info("Processing %i riboSeed runs; this can take a while",
                    len(riboSeed_jobs))

    pool = multiprocessing.Pool(processes=args.njobs)
    logger.debug("running the following commands:")
    logger.debug("\n".join(ribo_cmds))
    riboSeed_pool_results = [
        pool.apply_async(run_riboSeed_catch_errors,
                         (cmd,),
                         {"args": args,
                          "acc": acc,
                          "status_file": sfile,
                          "riboSeed_jobs": riboSeed_jobs})
        for acc, cmd, contigs, sfile, tax_d, _ in riboSeed_jobs]
    pool.close()
    pool.join()
    ribo_results_sum = sum([r.get() for r in riboSeed_pool_results])
    logger.debug("Sum of return codes (should be 0): %i", ribo_results_sum)

    for v in riboSeed_jobs:
        try:
            check_riboSeed_outcome(status_file, v[2])
            update_status_file(v[3], message="RIBOSEED COMPLETE")
            write_pass_fail(args, status="PASS", stage=v[0], note="")
            all_assemblies.append([v[2], v[4]])
        except riboSeedUnsuccessfulError as e:
            if v[4] == 1:
                write_pass_fail(args, status="ERROR",
                                stage=v[0],
                                note="riboSeed Error")
                logger.error(e)
            else:
                update_status_file(v[3], message="RIBOSEED COMPLETE")
                write_pass_fail(args, status="FAIL",
                                stage=v[0],
                                note="riboSeed unsuccessful")
                logger.error(e)

    #######################################################################
    extract16soutput = os.path.join(
        args.output_dir,
        "{}_ribo16s.fasta".format(args.organism_name.replace(" ", "_")))
    out_summary = os.path.join(args.output_dir, "sequence_summary.tab")
    for outf in [extract16soutput, out_summary]:
        if os.path.exists(outf):
            os.remove(outf)
    logger.info("attempting to extract 16S sequences for re-assemblies")
    n_extracted_seqs = 0
    singleline = True
    min_length = 1000  # minimum length seqeunce to extract
    for assembly, tax_d in all_assemblies:
        sra = str(Path(assembly).parents[4].name)
        barr_gff = os.path.join(args.output_dir, sra, "barrnap.gff")
        try:
            run_barrnap(assembly, barr_gff, logger)
            this_extracted_seqs = extract_16s_from_assembly(
                assembly, barr_gff, sra, extract16soutput, out_summary, args,
                singleline, tax_d, min_length, logger)
            n_extracted_seqs = n_extracted_seqs + this_extracted_seqs
        except extracting16sError as e:
            logger.error(e)
            write_pass_fail(args, status="ERROR", stage=sra,
                            note="unknown error extracting 16S")
        except barrnapError as e:
            logger.error(e)
            write_pass_fail(args, status="ERROR", stage=sra,
                            note="Error running barrnap")

    ###########################################################################
    logger.info("Wrote out %i sequences", n_extracted_seqs)
    if len(n_errors) != 0:
        logger.warning("Errors during run:")
        for k, v in n_errors.items():
            logger.warning("   " + k + " errors: " + str(v))
    if n_extracted_seqs == 0:
        write_pass_fail(args, status="FAIL", stage="global",
                        note="No 16s sequences detected in re-assemblies")
        logger.warning("No 16s sequences recovered. exiting")
        sys.exit()
    write_pass_fail(args, status="PASS", stage="global", note="")
    sys.exit()


if __name__ == '__main__':
    main()
