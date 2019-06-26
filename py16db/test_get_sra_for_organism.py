from run_all import download_SRA, filter_srapure
import os
import shutil
from nose.tools.nontrivial import with_setup


def test_filter_srapure():
    test_result = filter_srapure(path="./py16db/test_data/test_sraFind.txt",  organism_name="Lactobacillus oryzae", strains=1, get_all=True)
    assert ["DRR021662"] == test_result

def setup_func():
    "set up test fixtures"
    if os.path.exists("./py16db/test_function/"):
        shutil.rmtree("./py16db/test_function/")
    pass

def teardown_func():
    "tear down test fixtures"
    shutil.rmtree("./py16db/test_function/")


@with_setup(setup_func, teardown_func)


def test_download_SRA():
    os.makedirs("./py16db/test_function")
    download_SRA(destination="./py16db/test_function", SRA="SRR8443698")
    
    assert os.path.isdir("./py16db/test_function/raw"), "not creating correct subdirectory"
    
