from get_sra_for_organism import download_SRA, filter_srapure
import os
import shutil

def test_filter_srapure():
    test_result = filter_srapure("", "")
    print(test_result)
    assert ["DRR021662"] == filter_srapure("test_srapure", "Lactobacillus oryzae")

def test_download_SRA():
    if os.path.isdir("./test_function"):
       shutil.rmtree("./test_function")
    os.makedirs("./test_function", exist_ok=True)
    download_SRA(destination="./test_function", SRA="SRR8443698")
    assert os.path.isdir("./test_function/raw"), "not creating correct subdirectory"
    
