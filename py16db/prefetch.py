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

import logging as logger
from . import __version__

from py16db.FocusDBData import FocusDBData, fasterqdumpError
from py16db.shared_methods import filter_SRA


def get_args():
    parser = argparse.ArgumentParser(
        description="" )
    parser.add_argument("--focusDB_data", dest="focusDB_data",
                        help="path to data storage area; default ~/.focusDB/",
                        default=None)
    parser.add_argument("-n", "--organism_name",
                        help="genus or genus species in quotes",
                        required=True)
    parser.add_argument("-o", "--outdir",
                        help="output directory. This should be only be used" +
                        " to override the default cache; see" +
                        "https://github.com/ncbi/sra-tools/wiki/" +
                        "03.-Quick-Toolkit-Configuration",
                        required=False)
    parser.add_argument("--sraFind_path", dest="sra_path",
                        help="path to sraFind file; default is in ~/.focusDB/",
                        required=False)
    parser.add_argument("-S", "--n_SRAs", help="max number of SRAs to be run",
                        type=int, required=False)
    parser.add_argument("--seed",
                        help="random seed for subsampling referencese",
                        type=int, default=12345)
    parser.add_argument("--get_all",
                        help="if a biosample is associated with " +
                        "multiple libraries, default behaviour is to " +
                        "download the first only.  Use --get_all to " +
                        "analyse each library",
                        action="store_true", required=False)
    parser.add_argument("--output_cmds",
                        help="write commands to file, dont execute",
                        type=str, required=False)
    parser.add_argument("-b", "--batch_size",
                        help="batches for aspera calls",
                        default=10,
                        type=int, required=False)
    # We hide these options cause this script doesn't download the organisms
    parser.add_argument("--genomes_dir",
                        help=argparse.SUPPRESS)
    parser.add_argument("--prokaryotes", dest="prokaryotes",
                        help=argparse.SUPPRESS,
                        default=None,
                        required=False)
    return(parser.parse_args())


def make_prefetch_cmd(args, sras):
    cmd = "prefetch " + " ".join(sras)
    return cmd


def main():
    args = get_args()
    if shutil.which("prefetch") is None:
        raise ValueError("sra-tools must be installed and in PATH")
    fDB = FocusDBData(
        dbdir=None,
        refdir=args.genomes_dir,
        sraFind_data=args.sra_path,
        prokaryotes=args.prokaryotes)
    fDB.check_genomes_dir(org=args.organism_name)
    fDB.fetch_sraFind_data(logger=logger)

    filtered_sras = filter_SRA(
        sraFind=fDB.sraFind_data,
        organism_name=args.organism_name,
        strains=args.n_SRAs,
        thisseed=args.seed,
        logger=logger,
        get_all=args.get_all)
    fDB.fetch_sraFind_data(logger=logger)
    cmds = []
    for i in range(math.ceil(len(filtered_sras)/args.batch_size)):
        lowlm = i * args.batch_size
        uplm = min((i * args.batch_size) + args.batch_size, len(filtered_sras))
        these_sras = filtered_sras[lowlm: uplm]
        cmds.append(make_prefetch_cmd(args, these_sras))
    if args.output_cmds:
        with open(args.output_cmds, "w") as outf:
            for i, cmd in enumerate(cmds):
                if i % 2 == 0:
                    outf.write(cmd.replace("prefetch", "prefetch -t https") + "\n")
                else:
                    outf.write(cmd + "\n")
    else:
        ncmds = len(cmds)
        for i, cmd in enumerate(cmds):
            if i % 5 == 0:
                print("fetching %i of %i" % (i,  ncmds))
            subprocess.run(cmd,
                           shell=sys.platform != "win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)



if __name__ == '__main__':
    main()
