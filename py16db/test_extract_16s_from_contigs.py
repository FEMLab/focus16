from run_all import extract_16s_from_contigs

def test_extract():
    input_contigs="./data/test_data/mutans/downsampled/downsampledreadsf.fastq"
    barr_out="../test/0/riboSeed/barrnap"
    output="../test/0/riboSeed/ribo16s/"
    test_result=extract_16s_from_contigs(input_contigs=input_contigs,
                                         barr_out=barr_out, output=output)
    assert []  == test_result



    

