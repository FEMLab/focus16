import os
import shutil
import logging as logger

from pathlib import Path
from .run_all import parse_kraken_report

here = os.path.dirname(__file__)


testfile = os.path.join(here, "test_data", "kraken2.report")


def test_parse_kraken():
    full = parse_kraken_report(kraken2_report=testfile)


    ex_full = {'C': (94.51, '1236', 'Gammaproteobacteria'), 'G': (66.46, '561', 'Escherichia'), 'F': (89.02, '543', 'Enterobacteriaceae'), 'D': (95.73, '2', 'Bacteria'), 'S': (63.41, '562', 'Escherichia coli'), 'P': (94.51, '1224', 'Proteobacteria'), 'O': (92.68, '91347', 'Enterobacterales')}
    for k in ex_full.keys():
        assert full[k] == ex_full[k]
