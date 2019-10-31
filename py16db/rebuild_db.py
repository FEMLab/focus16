#!/usr/bin/env python3

import argparse
import subprocess
import sys
import shutil
import os
import subprocess
import math
import gzip
import math

import logging
logging.basicConfig(level="INFO")
from . import __version__

from py16db.FocusDBData import FocusDBData, fasterqdumpError
from py16db.shared_methods import filter_sraFind


def get_args():
    parser = argparse.ArgumentParser(
        description="rebuild the data db from SRAs in the data dir. " +
        "This is useful if something gets corrupted" )
    parser.add_argument("--focusDB_data", dest="focusDB_data",
                        help="path to data storage area; default ~/.focusDB/",
                        default=None)
    parser.add_argument("-o", "--outdir",
                        help="output directory. This should be only be used" +
                        " to override the default cache; see" +
                        "https://github.com/ncbi/sra-tools/wiki/" +
                        "03.-Quick-Toolkit-Configuration",
                        required=False)
    parser.add_argument("--sraFind_path", dest="sra_path",
                        help="path to sraFind file; default is in ~/.focusDB/",
                        required=False)
    # We hide these options cause this script doesn't download the organisms
    parser.add_argument("--genomes_dir",
                        help=argparse.SUPPRESS)
    parser.add_argument("--prokaryotes", dest="prokaryotes",
                        help=argparse.SUPPRESS,
                        default=None,
                        required=False)
    return(parser.parse_args())


def main():
    args = get_args()
    logger = logging.getLogger(__name__)
    if shutil.which("prefetch") is None:
        raise ValueError("sra-tools must be installed and in PATH")
    fDB = FocusDBData(
        dbdir=args.focusDB_data,
        refdir=args.genomes_dir,
        sraFind_data=args.sra_path,
        prokaryotes=args.prokaryotes,
        setup=False)
    fDB.fetch_sraFind_data(logger=logger)
    fDB.rebuild_fresh_db(logger=logger)



if __name__ == '__main__':
    main()
