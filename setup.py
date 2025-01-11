from setuptools import setup, find_packages
import os

def read_requirements():
    with open(os.path.join(os.path.dirname(__file__), 'requirements.txt')) as f:
        return f.read().splitlines()

setup(
    name='fits_extractor',
    version='0.1.0',
    author="Nicolas Obrier",
    author_email="nicolas.obrier@etu.unistra.fr",
    description='Library for extracting FITS headers and performing spatial queries',
    # long_description=open(os.path.join(os.path.dirname(__file__), "README.md")).read(),
    # long_description_content_type="text/markdown",
    url="",
    packages=find_packages(),  
    python_requires=">=3.6",
    install_requires=read_requirements(),  
)
