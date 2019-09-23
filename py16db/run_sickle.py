
import os
import sys
import shutil
import subprocess


def run_sickle(fastq1, fastq2, output_dir, run):
    """ run sickle for read trimming, return paths to trimmed reads
    """
    if shutil.which("sickle") is None:
        raise ValueError("sickle not found in PATH!")
    new_fastq1 = os.path.join(output_dir, "fastq1_trimmed.fastq")
    new_fastq2 = os.path.join(output_dir, "fastq2_trimmed.fastq")
    new_fastqs = os.path.join(output_dir, "singles_from_trimming.fastq")
    if fastq2 is None:
        ##since illumina 1.8, quality scores returned to sanger.
        cmd = str("sickle se -f {fastq1} " +
                  "-t sanger -o {new_fastq1}").format(**locals())
        new_fastq2 = None
    else:
        cmd = str("sickle pe -f {fastq1} -r {fastq2} -t sanger " +
                  "-o {new_fastq1} -p {new_fastq2} " +
                  "-s {new_fastqs}").format(**locals())
    if run:
        os.makedirs(output_dir)
        try:
            subprocess.run(cmd,
                           shell=sys.platform !="win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        except:
            cmd = cmd.replace("sanger", "solexa")
            try:
                subprocess.run(cmd,
                               shell=sys.platform !="win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
            except:
                raise ValueError("Error executing sickle cmd!")

    return (new_fastq1, new_fastq2)
