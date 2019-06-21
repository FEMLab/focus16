from run_all import run_riboseed
import os
import shutil

def test_riboseed():
    downreadsf="./data/test_data/mutans/downsampled/downsampledreadsf.fastq"
    downreadsr="./data/test_data/mutans/downsampled/downsampledreadsr.fastq"
    output_dir="./data/test_data/mutans/riboSeed/"  
    sra="./data/test_data/mutans/NC_004350.2.fna"
    if os.path.isdir("./data/test_data/mutans/riboSeed/"):
        shutil.rmtree("./data/test_data/mutans/riboSeed/")
        
    run_riboseed(sra=sra, readsf=downreadsf,
                 readsr=downreadsr, cores="4", threads="1",
                 v="1", output=output_dir)
    
    
    assert os.path.exists("./data/test_data/mutans/riboSeed/seed/final_long_reads/")

    return()
