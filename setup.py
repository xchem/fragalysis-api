from setuptools import setup, find_packages


#with open('requirements.txt') as f:
 #   requirements = f.read().splitlines()


setup(
    name='fragalysis_api',
    version='0.0.1',
  #  install_requires=requirements,
    author='Fragment 5',
    packages=find_packages(include=('fragalysis_api', 'fragalysis_api.*')),
    )
