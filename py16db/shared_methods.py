import random


def filter_SRA(sraFind, organism_name, strains, get_all, thisseed, logger):
    """sraFind [github.com/nickp60/srafind], contains"""
    results = []
    with open(sraFind, "r") as infile:
        for line in infile:
            split_line = [x.replace('"', '').replace("'", "") for x
                          in line.strip().split("\t")]
            if split_line[11].startswith(organism_name):
                if split_line[8].startswith("ILLUMINA"):
                    results.append(split_line[17])
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
