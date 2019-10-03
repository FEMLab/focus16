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


setup(name='focusDB',
      description='Draft genome reassembly using riboSeed, for the construction ' +
      'of high resolution 16S databases',
      version=verstr,
      url='https://github.com/FEMLab/focusdb',
      author='Ben Nolan,Nick Waters',
      author_email='N.BEN1@nuigalway.ie, nickp60@gmail.com',
      license='MIT',
      keywords='bioinformatics, assembly, 16s, database',
      packages=['py16db'],
      # scripts=['py16db/run_all.py',
      #          'py16db/test_alignment.py',
      #          'py16db/test_ave_read_len.py',
      #          'py16db/test_best_ref.py',
      #          'py16db/test_check_rDNA_num.py',
      #          'py16db/test_coverage.py',
      #          'py16db/test_extract_16s_from_contigs.py',
      #          'py16db/test_get_sra_for_organism.py',
      #          'py16db/test_parse_status.py',
      #          'py16db/test_run_sickle.py',
      #          'py16db/test_sralist.py',
      #          'py16db/get_n_genomes.py',
      #          'py16db/fetch_sraFind_data.py'],
      install_requires = ['biopython',  'plentyofbugs'],
      zip_safe=False,
      # include_package_data=True,
      # package_data={
      #     'test_data': [
      #         #'ecoli/*',
      #         # 'test_reads1.fq',
      #         # 'test_reads1.fq',
      #         'ribo16',
      #         'status/*',
      #         'test_16s_multilineSHORT.fasta',
      #         'test_sraFind.txt',
      #         'test_sralist.txt',
      #     ],
      # },
      entry_points={
          'console_scripts': [
              'focusDB = py16db.run_all:main',
              'focusDB-combine-with-silva = py16db.combine_focusdb_and_silva:main',
              'focusDB-shannon-entropy = py16db.calculate_shannon_entropy:main',
              'focusDB-align-and-trim = py16db.align_and_trim_focusdb:main',
              'focusDB-prefetch = py16db.prefetch:main'
          ],
      }
)
