import random


def filter_SRA(sraFind, organism_name, strains, get_all, thisseed, logger):
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

    if strains != 0:
        results = results[0:strains]
        logger.debug('Found SRAs: %s', results)

    sras = []
    for result in results:
        these_sras = result.split(",")
        if get_all:
            for sra in these_sras:
                sras.append(sra)
        else:
            sras.append(these_sras[0])

    return(sras)
