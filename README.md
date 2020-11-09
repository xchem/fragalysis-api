[![Build Status](https://travis-ci.org/xchem/fragalysis-api.svg?branch=master)](https://travis-ci.org/xchem/fragalysis-api)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/xchem/fragalysis-api.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/xchem/fragalysis-api/context:python)

Documentation: https://xchem.github.io/fragalysis-api/

# fragalysis-api

This api aims to allow any user to upload pdb files from the pdb or that they have created themselves, 
and analyse the ligand binding using the fragalysis webpage (https://fragalysis.diamond.ac.uk).

## Installation

Starting out by initialising an environment and activating it. 
Clone the repository and cd to the relevant directory. 
Install rdkit via conda, and the other dependencies via the setup.py file:
```bash
conda create -n fragalysis_env anaconda -y
conda activate fragalysis_env
conda install -c conda-forge rdkit -y

# Install our-bespoke version of gemmi # Required
git clone https://github.com/xchem/gemmi_pandda.git
cd gemmi_pannda/
pip install -U --force-reinstall .
cd ..

# Also Required
git clone https://github.com/xchem/pandda_gemmi.git
cd pandda_gemmi/
pip install -e .
cd ..

# Finally install the api
git clone "https://github.com/xchem/fragalysis-api.git"
cd fragalysis-api/
pip install -e .
cd ..
```

You can check if it has installed using:    ```conda list```

### How to use API

1. Set up environment
2. Download PDB files and query the PDB for structures of the same protein bound to the same or different ligands
3. Submit PDB files - you will be given a query ID 
4. Push your files into fragalysis and view them online :construction:
5. Analyse the binding of ligands to your target protein! :construction:

Other functionalities that are available:
- Import pdb files directly through the API
- Query the pdb for similar structures that also have ligands bound and have the option to import these structures
- Align files as individual monomers
- Align volume density files


### Enforced rules :scroll:

* The pdb file shall not be greater than 5mb.
* The pdb filename shall not contain non English language ascii characters 
and shall be between 4 and 20 characters in length.
* Each pdb file for alignment shall contain the same number of chains.
* All pdb files to be aligned must be in the same directory.
* If manually selecting a file for reference it must be in the same director as all pdb files for alignment. 
* PDB file must abide by best practices set out at https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html#note5
* In PDB file, ligands must be referred to by same code in 'HET' lines in header and in 'HETATM' lines in main part of file
* If providing cpp4 map files: All maps should be of the '.map' or '.cpp4' variety and labelled as:
    > ```
    > Aprot-x0001.pdb
    > Aprot-x0001_smiles.txt
    > Aprot_fofc.map
    > Aprot_2fofc.map
    > Aprot_event_0.map
* If you have multiple event maps increment the number and add accordingly.
* Maps shall be cut/masked by yourself using your tool of choice after alignment (sorry)

### 2.1 How to download PDB files

Move to the fragalysis-api/xcimporter directory. 

You will need two bits of information:
1. Your 'user ID' - this is your name followed by your protein name. For example, Anna looking at protein ATAD2 wil have username 'Anna_ATAD2'.
2. The PDB code you would like to download. For example, '6epu'.

To download the PDB file 
```
python pdbimporter.py -id [user_id] -pdb [pdb code]
```
In our example, this would be 
```python pdbimporter.py -id Anna_ATAD2 -pdb 6epu```

Alternatively, you can upload your own PDB files. :construction:

### 2.2 How to query the PDB file

Note: you don't need to manually download a pdb file before querying the PDB for structures of the same protein

Similar to downloading PDB files, however initially all you need is the PDB code and the chain you would like to query. The command is 
```
python pdbquery.py -pdb [pdb code] -chain [required chain]
```
The API will then query the PDB for structures of the same protein that also have ligands bound. You will be asked if you would like to see a list of these structures and ligands, and if you would like to download all of these pdb files in bulk. If you choose to download these pdb files, you will be asked for your user id. This is as before: your name followed by the protein name (e.g. Anna_ATAD2).

### 3.How to submit PDB files for conversion to a fragalysis friendly format (fff)

Once you have the files downloaded, they need to be processed before they can be visualised in fragalysis. This is done using
```
python xcimporter -id [user id] 
```
Default directories are used. These can however be changed by using ```-i [input directory] ``` or   ```-o [output directory] ``` if this is required.

For example a more elaborate conversion would be:
```
pythom xcimporter -i [input directory] -o [output directory] -t [target name] -m 
```
Which will align the input files with respects to individual chains inside the .pdb files (`-m`) and save the output
to a folder specified by `-t [target name]`.

The terminal will let you know when the conversion has been successful and if there are any files that have been found to be incompatible with the API. We are working to minimize any incompatibilities.

The expected output of the xcimporter is a folder located at `[outputdirectory]/[target name]` which will
contain two main folders `aligned` and `crystallographic`. The `aligned` folder will contain multiple sub-folders, one per ligand. These will be labelled as the name of the `pdbfile_[lig_number]` or `pdbfile_[lig_number][chain_letter]` if monomer mode was specified.
Inside each of these subfolders contains all the information that is required by the fragalysis platform to display and elaborate a given ligand.
The `crystallographic` folder contains an exact copy of the data supplied to the above command.

### Who we are
We are the Fragment 5, a group of students at the University of Oxford.

Anna :whale:   
Maranga :fire:  
George a.k.a Joff Boff :poop:  
Tobias :cow:  
Alister :panda_face:  

We are looked after by Rachael :crown:. 

Finally, this project has been ruined by Tyler.




