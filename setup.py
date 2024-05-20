from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from rnavirhost import __version__

setup(name='rnavirhost',
      version=__version__,
      packages=find_packages(),
      package_data={"rnavirhost":["virus/*", "model/*"]},
      #install_requires=['blast>=1.0.1'],
      description='RNAVirHost: a machine learning-based method for predicting hosts of RNA viruses through viral genomes',
      keywords=['Bioconda RNAvirus HostPrediction MachineLearning Metagenomics'],
      classifiers=[],
      url='https://github.com/GreyGuoweiChen/VirHost.git',
      author='CHEN Guowei',
      author_email='gwchen3-c@my.cityu.edu.hk',
      license='MIT',
      install_requires=[
#        'blast>=2.12.0',
        'numpy>=1.23.5',
        'pandas>=2.0.3',
        'scikit-learn>=1.1.3',
        'xgboost>=1.7.4',
        'biopython>=1.83',
#        'prodigal>=2.6.3'
    ],
      entry_points={
        "console_scripts":[
        "rnavirhost = rnavirhost.rnavirhost:main"
#        "rnavirhost = rnavirhost.rnavirhost:main",
#        "determine_order = rnavirhost.determine_order:main",
        ]},
      include_package_data=True,
      zip_safe=False
)
