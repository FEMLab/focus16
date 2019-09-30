import os
from .run_all import parse_status_file, update_status_file
here = os.path.dirname(__file__)
fpath = os.path.join(here, "sample_status_file")

with open(fpath, "w") as f:
    f.write("FIRST THING")

def test_parse_status_file():
    statuses = parse_status_file(fpath)
    assert statuses == ["FIRST THING"]


def test_update_status_file():
    update_status_file(fpath, message="DONE")
    statuses = parse_status_file(fpath)
    update_status_file(fpath, to_remove = ["DONE"], message="THIRD THING")
    statuses = parse_status_file(fpath)
    print(statuses)
    assert statuses.sort() == ["FIRST THING", "THIRD THING"].sort()
