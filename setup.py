from setuptools import setup
import re

with open("README.md", "r") as fh:
    long_description = fh.read()
    
VERSIONFILE = "py16db/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))


setup(name='16db',
      description='Draft genome reassembly using riboSeed, for the construction ' + 
      'of high resolution 16S databases',
      version=verstr,
      url='http://github.com/BenNolann/16db',
      author='Ben Nolan',
      author_email='N.BEN1@nuigalway.ie',
      license='MIT',
      keywords='bioinformatics, assembly, 16s, database',
      packages=['py16db'],
      scripts=['py16db/run_all.py',
               'py16db/test_ave_read_len.py',
               'py16db/test_best_ref.py',
               'py16db/test_coverage.py',
               'py16db/test_get_sra_for_organism.py',
               'py16db/test_run_sickle.py',
               'py16db/get_n_genomes.py',
               'py16db/fetch_sraFind_data.py'],
      install_requires = ['biopython'],
      zip_safe=False,
      include_package_data=True,
      package_data={
        'example_data':
            ['example_data/NC_013928.1.fna',
             'example_data/reads1.fq.gz',
             'example_data/reads2.fq.gz',
             'example_data/README.md',
             ],
        'test_data':
            ['plasmids/*',
             'test_reads1.fq',
             'sraFind-Contig-biosample-with-SRA-hits.txt',
             'NC_017659.1.fasta',
             'test_16s_multiline.fasta',
             ],
        },
      entry_points={ 
        'console_scripts': [
            '16db = py16db.run_all:main',
            ],
        }
      )
                    
          
