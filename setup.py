from setuptools import setup, find_packages
import os


def read_requirements():
    with open('requirements.txt') as f:
        return f.read().splitlines()

setup(
    name='fits_extractor',
    version='0.1',
    description='Library for extracting FITS headers and performing spatial queries',
    packages=find_packages(),  
    install_requires=read_requirements(),  
)
