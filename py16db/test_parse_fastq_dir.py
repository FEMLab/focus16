import os
from pathlib import Path
import logging as logger
from .run_all import check_fastq_dir
here = os.path.dirname(__file__)

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
for k, v in thesedirs.items():
    base = os.path.join(testdir, k)
    os.makedirs(base, exist_ok=True)
    for d in v[0]:
        Path(os.path.join(base, "test" + d)).touch()


def test_parse_status_file():
    for k, v in thesedirs.items():
        base = os.path.join(testdir,  k)
        if k.startswith("bad"):
            fwd, rev, mgs = check_fastq_dir(base, logger)
            assert mgs != "", "no error thrown"
        else:
            fwd, rev, mgs = check_fastq_dir(base, logger)
            if v[1]:
                assert fwd == os.path.join(base, "test_1.fastq"), "missing forward library"
            if v[2]:
                assert rev == os.path.join(base, "test_2.fastq"), "missing rev library"
            if v[3]:
                assert fwd == os.path.join(base, "test.fastq"), "missing single library"
