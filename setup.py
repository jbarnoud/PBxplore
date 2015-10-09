#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, absolute_import

from setuptools import setup, find_packages
import os

import pbxplore

here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(here, 'README.md')) as f:
    readme = f.read()


setup(

    name='pbxplore',
    version=pbxplore.__version__,

    description="PBxplore is a suite of tools dedicated to Protein Block analysis.",
    long_description=readme,

    url='https://github.com/pierrepo/PBxplore',

    # Author details
    author='Pierre Poulain',
    author_email='pierre.poulain@cupnet.net',

    license='MIT',

    classifiers=[
        'Development Status :: 4 - Beta',

        # Indicate who your project is intended for
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Physics',

        'License :: OSI Approved :: MIT License',

        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
    ],

    install_requires=['numpy'],
    packages=find_packages(exclude=['doc']),
    entry_points={
        'console_scripts': [
            'PBassign = pbxplore.scripts.PBassign:pbassign_cli',
            'PBclust  = pbxplore.scripts.PBclust:pbclust_cli',
            'PBcount  = pbxplore.scripts.PBcount:pbcount_cli',
            'PBstat   = pbxplore.scripts.PBstat:pbstat_cli',
        ],
    },

)