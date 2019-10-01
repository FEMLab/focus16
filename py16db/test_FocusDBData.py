import os
import shutil
from pathlib import Path
import logging as logger
from .FocusDBData import FocusDBData
here = os.path.dirname(__file__)
thisdbdir = os.path.join(here, ".test_db_dir", "")


if os.path.exists(thisdbdir):
    shutil.rmtree(thisdbdir)

def test_make_data_dir():
    fDB = FocusDBData(dbdir = thisdbdir)
    assert os.path.exists(thisdbdir), "error creating dir"

def test_start_SRAs():
    fDB = FocusDBData(dbdir = thisdbdir)
    assert len(fDB.SRAs.keys()) == 0, "list not starting empty to begin with!"

def test_add_SRAs():
    fDB = FocusDBData(dbdir = thisdbdir)
    fDB.update_manifest(newacc="SRRtest123", newstatus="mostly_ok",
                        organism = "test", logger=logger)
    print(fDB.SRAs)
    assert fDB.SRAs['SRRtest123']['status'] == "mostly_ok", \
        "error adding SRA to list"
