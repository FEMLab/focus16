import random
import os
import sys
import subprocess
from Bio import SeqIO
import gzip


def filter_sraFind(sraFind, organism_name, strains, get_all, thisseed,
               use_available, logger):
    results = get_lines_from_sraFind(sraFind, organism_name)
    random.seed(thisseed)
    random.shuffle(results)
    # log  a sane amount
    if len(results) < 20:
        logger.debug('All matching SRAs from sraFind: %s', results)
    else:
        logger.debug('All matching SRAs from sraFind: %s and %i others', results[0:20], len(results)-20)

    if use_available:
        focusDB_dir = os.path.dirname(sraFind)
        # if dir exists
        possible_sras  = [x for x in results if os.path.exists(os.path.join(focusDB_dir, x))]
        # if data in dir
        possible_sras  = [x for x in  possible_sras if os.listdir(os.path.join(focusDB_dir, x))]
        logger.debug("using locally available SRAs; found %i", len(possible_sras))
        if strains == 0:
            pass  # just use  as many available as possible
        else:
            # in case we dont have enough locally, make an "set" from (the ones
            # we have extended with the ones we dont), and get the first N;
            # this puts tthe local ones first in the list.
            # dont rrefactor with a set -- they are unordered
            possible_sras.extend([x for x in results if x not in possible_sras])
            possible_sras = possible_sras[0: strains]
        logger.debug(possible_sras)
    else:
        if strains != 0:
            possible_sras = results[0:strains]
        else:
            possible_sras = results
    logger.debug('Selected the following SRAs: %s', possible_sras)

    sras = []
    for result in possible_sras:
        these_sras = result.split(",")
        if get_all:
            for sra in these_sras:
                sras.append(sra)
        else:
            sras.append(these_sras[0])
    logger.info('Processing the following SRAs: %s', sras)
    return(sras)


def get_lines_from_sraFind(sraFind, organism_name):
    """sraFind [github.com/nickp60/srafind], contains"""
    results = []
    with open(sraFind, "r") as infile:
        for i, line in enumerate(infile):
            if i == 0:
                header = line.strip().replace('"', '').replace("'", "").split("\t")
                org_col = header.index("organism_ScientificName")
                plat_col  = header.index("platform")
                run_SRAs = header.index("run_SRAs")
            split_line = [x.replace('"', '').replace("'", "") for x
                          in line.strip().split("\t")]
            if split_line[org_col].startswith(organism_name):
                if split_line[plat_col].startswith("ILLUMINA"):
                    results.append(split_line[run_SRAs])
    results = [x for x in results if x != ""]
    return results


def get_ave_read_len_from_fastq(fastq1, logger=None):
    """return average read length in fastq1 file from first N reads
    """
    nreads = None
    nreads = 1000
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
    return float(tot / i)


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
    logger.debug("processing %s", assembly)
    tax_string = tax_d["S"][2]
    # if no label at some level, give the next one
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
            # need this: catches index errors from the comment lines
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
                if (end - start) < min_length:
                    logger.info("ignoring short rRNA:")
                    logger.info(line)
                    continue
                thisid = "{}_{}".format(sra, rrn_num)
                results16s[thisid] = [chrom, start, end, line[6]]
                with open(assembly, "r") as asmb:
                    for rec in SeqIO.parse(asmb, "fasta"):
                        #logger.debug("%s\t%s", chrom, rec.id)
                        if rec.id != chrom:
                            continue
                        seq = rec.seq[start + 1: end + 1]
                        if ori == "-":
                            seq = seq.reverse_complement()
                        thisidcoords = "{thisid}.{start}.{end}".format(
                            **locals())
                        # Need to disable linewrapping for use with SILVA, etc
                        if singleline:
                            seqstr = str(seq.transcribe())
                            outf.write(
                                str(">{thisidcoords} {big_tax_string}\n" +
                                    "{seqstr}\n").format(**locals()))
                        else:
                            SeqIO.write(
                                SeqRecord(
                                    seq.transcribe(), id=thisidcoords,
                                    description=big_tax_string),
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
