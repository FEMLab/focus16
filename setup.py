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
      install_requires = ['biopython',  'plentyofbugs'],
      zip_safe=False,
      entry_points={
          'console_scripts': [
              'focusDB = py16db.run_focusDB:main',
              'focusDB-combine-with-silva = py16db.combine_focusdb_and_silva:main',
              'focusDB-shannon-entropy = py16db.calculate_shannon_entropy:main',
              'focusDB-align-and-trim = py16db.align_and_trim_focusdb:main',
              'focusDB-prefetch = py16db.prefetch:main',
              'focusDB-rebuild = py16db.rebuild_db:main'
          ],
      }
)
