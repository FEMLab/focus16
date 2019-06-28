from .run_sickle import run_sickle
import os
import shutil
from Bio import SeqIO
from nose.tools.nontrivial import with_setup


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def setup_func():
    "set up test fixtures"
    if os.path.exists("./py16db/testsickle/"):
        shutil.rmtree("./py16db/testsickle/")
    pass

def teardown_func():
    "tear down test fixtures"
    shutil.rmtree("./py16db/testsickle/")


@with_setup(setup_func, teardown_func)
def test_sickle_singles():
    sickle_test_dir = "./py16db/testsickle/"
    fastq1 = os.path.join(os.path.dirname(__file__), "test_data", "mutans", "downsampled", "downsampledreadsf.fastq")
    fastq2 = None
    if os.path.exists(sickle_test_dir):
        shutil.rmtree(sickle_test_dir)
    new_fastq1, new_fastq2 = run_sickle(fastq1=fastq1, fastq2 = None,
                                        output_dir=sickle_test_dir)
    assert new_fastq2 == None
    assert file_len(fastq1) > file_len(new_fastq1)
