#!/usr/bin/env python

try:
    from setuptools import setup
except ImportError:
    from ez_setup import use_setuptools
    use_setuptools()
    from setuptools import setup

import os


__version__ = '0.1'

setup(
    name='jwstpipett',
    version=__version__,
    description='James Webb Pipeline Testing Tools',
    author='Patrick Kavanagh',
    author_email='pkavanagh@cp.dias.ie',
    install_requires=[
        #'configobj',
        #'jwst',
        #'logging',
    ],
    keywords = ['testing'], 
    packages=['jwstpipett'],
    scripts=["jwstpipett/plot_ramps.py"],
    #package_data={'mirie2e': ['cfg_files/*.cfg']},
    classifiers = [
        "Programming Language :: Python :: 3.6",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research/DataAnalysis",
        #"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Topic :: Scientific/Engineering"
    ]#,
    #include_package_data=True
)
