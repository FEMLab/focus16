#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import sys
import os
import subprocess
import shutil
import logging
import glob
import multiprocessing
import tempfile


from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


from py16db.shared_methods import filter_sraFind, \
    extract_16s_from_assembly, run_barrnap, parse_kraken_report, run_kraken2




def setup_logging(args):  # pragma: nocover
    if (args.verbosity * 10) not in range(10, 60, 10):
        raise ValueError('Invalid log level: %s' % args.verbosity)
    logging.basicConfig(
        level=logging.DEBUG,
        filemode='w',
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
        description="this utilty annotates taxonomy from assembly" +
        " and extracts the 16s, outputting to match silva's format." +
        " Not part of the focusDB pipeline, but used in manuscript to assess" +
        "SPAdes de novo assembly.",
        add_help=False)  # to allow for custom help
    parser.add_argument(
        "contigs",
        help="input contigs")
    parser.add_argument(
        "outfile",
        help="output file")
    parser.add_argument(
        "--name",
        help="sequence naming prefix",
        default="focusDB")
    parser.add_argument(
        "--focusDB_data", dest="focusDB_data",
        help="path to data storage area; default ~/.focusDB/",
        default=None)
    parser.add_argument(
        "--kraken2_dir", action="store",
        help="path to kraken dir; default is " +
        "in ~/.focusDB/.",
        default=os.path.expanduser("~/.focusDB/"),
        required=False)
    parser.add_argument(
        "--clobber", action="store_true",
        help="destroy existing output")
    parser.add_argument(
        "--kraken_mem_mapping", action="store_true",
        help="use this flag to load kraken2 db via disk " +
        "instead of RAM for taxonomic assignment. " +
        "automatically enabled if --memory < 20GB",
        required=False)
    parser.add_argument(
        "--njobs",
        help="how many jobs to run concurrently " +
        "via multiprocessing or SGE. --cores and --memory is per job",
        default=1, type=int)
    parser.add_argument(
        "--memory",
        help="PER JOB: amount of RAM to be used. riboSeed " +
        "needs 10GB ram to run optimally; less, riboSeed " +
        "runs in --serialize mode to prevent memory errors" +
        "during subassemblies",
        default=4,
        required=False, type=int)
    parser.add_argument(
        "--threads",
        action="store",
        default=1, type=int,
        choices=[1, 2, 4],
        help="if your cores are hyperthreaded, set number" +
        " threads to the number of threads per processer." +
        "If unsure, see 'cat /proc/cpuinfo' under 'cpu " +
        "cores', or 'lscpu' under 'Thread(s) per core'." +
        ": %(default)s")
    parser.add_argument(
        "--cores",
        action="store",
        default=1, type=int,
        help="number of cores to run kraken with" +
        ": %(default)s")
    parser.add_argument(
        "-v", "--verbosity", dest='verbosity',
        action="store",
        default=2, type=int, choices=[1, 2, 3, 4, 5],
        help="Logger writes debug to file in output dir; " +
        "this sets verbosity level sent to stderr. " +
        " 1 = debug(), 2 = info(), 3 = warning(), " +
        "4 = error() and 5 = critical(); " +
        "default: %(default)s")
    args = parser.parse_args()
    return(args)


def main():
    args = get_args()
    logger = setup_logging(args)
    tmp = tempfile.mkdtemp()
    logger.debug("Logging to %s", tmp)
    kraken2_report_output = os.path.join(tmp, "kraken2.report")
    #######################################################################
    report_output = run_kraken2(
        args,
        contigs=args.contigs,
        dest_prefix=os.path.join(tmp, "kraken2"),
        db=os.path.join(args.kraken2_dir, "minikraken2_v2_8GB_201904_UPDATE", ""),
        logger=logger)
    tax_dict = parse_kraken_report(kraken2_report=report_output)
    print(report_output)
    barr_gff = os.path.join(tmp, "barrnap_full.gff")
    subprocess.run("barrnap %s > %s" %(args.contigs, barr_gff),
                   shell=sys.platform !="win32",
                   stdout=subprocess.PIPE,
                   stderr=subprocess.PIPE,
                   check=True)
    this_extracted_seqs = extract_16s_from_assembly(
        args.contigs, barr_gff, sra=args.name,
        output=args.outfile, output_summary=os.path.join(tmp, "summary"), args=args,
        singleline=True, tax_d=tax_dict, min_length=1200, logger=logger)

    sys.exit()


if __name__ == '__main__':
    main()
