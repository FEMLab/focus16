from setuptools import setup

setup(name='16db',
      version='0.1',
      description='Draft genome reassembly using riboSeed, for the construction ' + 
      'of new databases',
      url='http://github.com/BenNolann/16db',
      author='Ben Nolan',
      author_email='N.BEN1@nuigalway.ie',
      license='MIT',
      packages=['py16db'],
      scripts=['py16db/get_sra_for_organism.py', 'best_ref.py', 
               'get_n_genomes.py', 'riboseed.py'],
      install_requires = ["biopython"],
      zip_safe=False)
