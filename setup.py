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
      url='https://github.com/GreyGuoweiChen/VirHost.git',
      author='CHEN Guowei',
      author_email='gwchen3-c@my.cityu.edu.hk',
      license='MIT',
      entry_points={
        "console_scripts":[
        "RNAVirHost = rnavirhost.RNAVirHost:main"
#        "rnavirhost = rnavirhost.rnavirhost:main",
#        "determine_order = rnavirhost.determine_order:main",
        ]},
      include_package_data=True,
      keywords=[],
      zip_safe=False)
