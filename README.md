[![Build Status](https://travis-ci.org/xchem/fragalysis-api.svg?branch=master)](https://travis-ci.org/xchem/fragalysis-api)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/xchem/fragalysis-api.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/xchem/fragalysis-api/context:python)

# fragalysis-api

This api aims to allow any user to upload pdb files from the pdb or that they have created themselves, and analyse the ligand binding using the fragalysis webpage (https://fragalysis.diamond.ac.uk).

### How to use API

(note: not live yet)

1. Set up environment
2. Submit PDB files - you will be given a query ID 
3. Push your files into fragalysis and view them online
4. Analyse the binding of ligands to your target protein!

Other functionalities that are available:
- Import pdb files directly through the API
- Query the pdb for similar structures that also have ligands bound and have the option to import these structures

##### How to set install, update and save the fragalysis_env conda environment.

This command exports the env you are in to a .yml file.
```
conda env export > environment.yml
```
This command creates the env from the .yml file.
```
conda env create -f environment.yml
```
This command updates the current env you are in with the .yml file.
```
conda env update -f environment.yml
```
This command activates the created environment 
```
conda activate fragalysis_env
```
To install setup.py make sure you  are in the fragalysis_env: 
```
pip install -e .
```
You can check if it has installed using
```
conda list
```


### Enforced rules :scroll:
* The pdb file shall not be greater than 5mb.
* The pdb filename shall not contain non English language ascii characters and shall be between 4 and 20 characters in length.
* Each pdb file for alignment shall contain the same number of chains.
* All pdb files to be aligned must be in the same directory.
* If manually selecting a file for reference it must be in the same director as all pdb files for alignment. 

### Who we are

We are the Fragment 5, a group of students at the University of Oxford.

Anna :whale:   
Maranga :fire:  
George a.k.a Joff Boff :squirrel:  
Tobias :cow:  
Alister :panda_face:  

We are looked after by Rachael :crown:

