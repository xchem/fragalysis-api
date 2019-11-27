from setuptools import setup, find_packages

setup(
    name='fragalysis_api',
    version='0.0.1',
    author='Fragment 5',
    packages=find_packages(include=('fragalysis_api', 'fragalysis_api.*')),
    )
