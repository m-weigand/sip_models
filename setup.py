#!/usr/bin/env python
from setuptools import setup
version_short = '0.1'
version_long = '0.1.3'

extra = {}

if __name__ == '__main__':
    setup(
            name='sip_models',
            version=version_long,
            description=''.join((
                'Spectral induced polarization (SIP) models based on ',
                'the Cole-Cole model',
            )),
            author='Maximilian Weigand',
            author_email='mweigand@geo.uni-bonn.de',
            url='https://github.com/m-weigand/sip_models',
            package_dir={'': 'lib'},
            packages=[
                'sip_models',
                'sip_models.res',
                'sip_models.cond',
            ],
            install_requires=[
                'numpy',
                'scipy>=0.12',
                'matplotlib'
            ],
            classifiers=[
                "Development Status :: 4 - Beta",
                "License :: OSI Approved :: GNU Lesser General " +
                "Public License v3 (LGPLv3)",
                "Programming Language :: Python :: 3.4",
                "Programming Language :: Python :: 3.5",
                "Programming Language :: Python :: 3.9",
                "Intended Audience :: Science/Research",
            ],
            **extra
    )
