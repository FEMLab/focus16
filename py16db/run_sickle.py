import os
import sys
import shutil
import subprocess

def run_sickle(fastq1, fastq2, output_dir):
    """ run sickle for read trimming, return paths to trimmed reads
    """
    if shutil.which("sickle") is None:
        raise ValueError("sickle not found in PATH!")
    os.makedirs(output_dir)
    new_fastq1 = os.path.join(output_dir, "fastq1_trimmed.fastq")
    new_fastq2 = os.path.join(output_dir, "fastq2_trimmed.fastq")
    new_fastqs = os.path.join(output_dir, "singles_from_trimming.fastq")
    if fastq2 is None:
        cmd = str("sickle se -f {fastq1} " +
                  "-t illumina -o {new_fastq1}").format(**locals())
        new_fastq2 = None
    else:
        cmd = str("sickle pe -f {fastq1} -r {fastq2} -t illumina " +
                  "-o {new_fastq1} -p {new_fastq2} " +
                  "-s new_fastqs}" +
                  "sickle pe -f {fastq1}").format(**locals())
    print(cmd)
    try:
        subprocess.run(cmd,
                       shell=sys.platform !="win32",
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE,
                       check=True)
    except:
        cmd = cmd.replace("illumina", "solexa")
        try:
            subprocess.run(cmd,
                           shell=sys.platform !="win32",
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True)
        except:
            cmd = cmd.replace("solexa", "sanger")
            try:
                subprocess.run(cmd,
                               shell=sys.platform !="win32",
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               check=True)
            except:
                print(cmd)
                raise ValueError("Error executing sickle cmd!")

    return (new_fastq1, new_fastq2)
