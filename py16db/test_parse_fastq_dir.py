import os
import shutil
from pathlib import Path
import logging as logger

here = os.path.dirname(__file__)
from py16db.FocusDBData import FocusDBData, fasterqdumpError

thesedirs = {
    # dir,      REads,                  should have F,   R, S,  message
    # note we reject single
    "goodA": [["_1.fastq", ".fastq", "_2.fastq"],
              True, True, False, ""],
    "goodB": [["_1.fastq", "_2.fastq"],
              True, True, False,  ""],
    "goodC": [[".fastq"],
              False, False, True, ""],
    "badA":  [["_2.fastq"],
              False, False, False, "error"],
    "badB":  [[".fastq", "_2.fastq"],
              False, False, False, "f"]
}


testdir = os.path.join(here, "tmp_fastqs")
if os.path.exists(testdir):
    shutil.rmtree(testdir)
for k, v in thesedirs.items():
    base = os.path.join(testdir, k)
    os.makedirs(base, exist_ok=True)
    for d in v[0]:
        Path(os.path.join(base, "test" + d)).touch()


def test_fastq_files():
    fDB = FocusDBData(dbdir = testdir,
                      refdir = testdir + "res",
                      prokaryotes="proks",
                      sraFind_data="sraFind.txt")

    for k, v in thesedirs.items():
        base = os.path.join(testdir,  k)
        if k.startswith("bad"):
            fwd, rev, mgs = fDB.check_fastq_dir(base, logger)
            assert mgs != "", "no error thrown"
        else:
            fwd, rev, mgs = fDB.check_fastq_dir(base, logger)
            if v[1]:
                assert fwd == os.path.join(base, "test_1.fastq"), "missing forward library"
            if v[2]:
                assert rev == os.path.join(base, "test_2.fastq"), "missing rev library"
            if v[3]:
                assert fwd == os.path.join(base, "test.fastq"), "missing single library"
