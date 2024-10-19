from setuptools import setup, find_packages
import os

# Lire le fichier requirements.txt
def read_requirements():
    with open('requirements.txt') as f:
        return f.read().splitlines()

setup(
    name='fits_extractor',
    version='0.1',
    description='Library for extracting FITS headers and performing spatial queries',
    packages=find_packages(),  # Inclure les paquets automatiquement
    install_requires=read_requirements(),  # Utiliser requirements.txt pour les dépendances
)
