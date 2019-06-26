from setuptools import setup

setup(name='16db',
      version='0.1.0',
      description='Draft genome reassembly using riboSeed, for the construction ' + 
      'of new 16S databases',
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
               'py16db/test.py',
               'get_n_genomes.py'],
      install_requires = ["biopython"],
      zip_safe=False,
      include_package_data=True,
      package_data={
          'py16db':[
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
                   ],
              },
          entry_points={ 
              'console_scripts': [
                  '16db = py16db.run_all.py:main',
                  ],
              }
          )
                    
          
