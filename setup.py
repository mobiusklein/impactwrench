from distutils.core import setup
from setuptools import setup

def readme():
    with open('README.rst') as f:
         return f.read()

setup(name = 'impactwrench',
      version = '0.1',
      description ='ImpactWrench is a Python CLI application for LC-MS and LC-MS/MS data analysis  ',
      author = 'Kundai Sachikonye',
      author_email = 'k.sachikonye@uke.de',
      scripts = ['bin/impactwrench'],
      license = 'MIT',
      install_requires =['pymzml', 'scipy', 'numpy','matplotlib','ursgal', 'pymc3', 'theano','seaborn','sklearn', 'mlxtend'],
      test_suite = 'nose.collector',
      tests_require = ['nose', 'nose-cover3'],
      include_package_data = True,
      zip_safe= False )