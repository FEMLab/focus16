import random
import os

def filter_SRA(sraFind, organism_name, strains, get_all, thisseed,
               use_available, logger):
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
    random.seed(thisseed)
    random.shuffle(results)
    # log  a sane amount
    if len(results) < 20:
        logger.debug('All matching SRAs from sraFind: %s', results)
    else:
        logger.debug('All matching SRAs from sraFind: %s and %i others', results[0:20], len(results)-20)

    if use_available:
        focusDB_dir = os.path.dirname(sraFind)
        possible_sras  = [x for x in results if os.path.exists(os.path.join(focusDB_dir, x))]
        logger.debug("using locally available SRAs; found %i", len(possible_sras))
    else:
        possible_sras = results
    if strains != 0:
        possible_sras = possible_sras[0:strains]
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
    return([x for x in sras if  x != ""])
