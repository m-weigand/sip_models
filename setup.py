#!/usr/bin/env python
from setuptools import setup
import subprocess
import os
import sys
# from setuptools import find_packages
# find_packages

# under windows, run
# python.exe setup.py bdist --format msi
# to create a windows installer

version_short = '0.1'
version_long = '0.1.0'
# if we are in a git directory, use the last git commit as the version
cmd = 'git log -1 --format=%H'
try:
    if os.path.isdir('.git'):
        git_output = subprocess.check_output(cmd, shell=True).strip()
        version_long += '+' + git_output
except:
    pass

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
          url='http://www.geo.uni-bonn.de/~mweigand',
          # find_packages() somehow does not work under Win7 when creating a
          # msi # installer
          # packages=find_packages(),
          package_dir={'': 'lib'},
          packages=['sip_models', 'sip_models.res' ],
          # scripts=['src/dd_single/dd_single.py',
          #          'src/dd_time/dd_time.py',
          #          'src/dd_space_time/dd_space_time.py',
          #          'src/dd_test/dd_test.py',
          #          'src/ddps/ddps.py',
          #          'src/ddpt/ddpt.py',
          #          'src/ddpst/ddpst.py',
          #          'src/ddplot/ddplot.py'],
          install_requires=['numpy', 'scipy>=0.12', 'matplotlib'],
          **extra
          )
