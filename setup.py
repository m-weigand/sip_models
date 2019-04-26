#!/usr/bin/env python
from setuptools import setup
import sys
# from setuptools import find_packages
# find_packages

# under windows, run
# python.exe setup.py bdist --format msi
# to create a windows installer

version_short = '0.1'
version_long = '0.1.2'

extra = {}
if sys.version_info >= (3,):
    print('V3')
    extra['use_2to3'] = True

if __name__ == '__main__':
    setup(name='sip_models',
          version=version_long,
          description='Spectral induced polarization (SIP) models based on ' +
          'the Cole-Cole model',
          author='Maximilian Weigand',
          author_email='mweigand@geo.uni-bonn.de',
          url='https://github.com/m-weigand/sip_models',
          # find_packages() somehow does not work under Win7 when creating a
          # msi # installer
          # packages=find_packages(),
          package_dir={'': 'lib'},
          packages=[
              'sip_models',
              'sip_models.res',
              'sip_models.cond',
          ],
          # scripts=[,],
          install_requires=['numpy', 'scipy>=0.12', 'matplotlib'],
          classifiers=[
              "Development Status :: 4 - Beta",
              "License :: OSI Approved :: GNU Lesser General " +
              "Public License v3 (LGPLv3)",
              "Programming Language :: Python :: 3.4",
              "Intended Audience :: Science/Research",
          ],
          **extra
          )
