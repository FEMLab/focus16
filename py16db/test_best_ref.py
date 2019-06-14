from run_all import pob
import os
import shutil

def test_pob():
    plasmids = "./plentyofbugs/test_data/plasmids/" 
    reads = "./plentyofbugs/test_data/test_reads1.fq" 
    output_dir="./pob_test_result/"
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    test_result = pob(genomes_dir=plasmids, readsf=reads, output_dir=output_dir)
    print(test_result)
    assert ["./plentyofbugs/test_data/plasmids/NC_009837.1.fasta", "0.0340249"] == test_result
    return()

    
