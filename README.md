# fragalysis-api

[![build main](https://github.com/InformaticsMatters/fragalysis-api/actions/workflows/build-main.yaml/badge.svg)](https://github.com/InformaticsMatters/fragalysis-api/actions/workflows/build-main.yaml)

[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/xchem/fragalysis-api.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/xchem/fragalysis-api/context:python)

Documentation: https://xchem.github.io/fragalysis-api/

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
cd gemmi_pandda/
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
* If manually selecting a file for reference it must be in the same directory as all pdb files for alignment. 
* PDB file must abide by best practices set out at https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html#note5
* In PDB file, ligands must be referred to by same code in 'HET' lines in header and in 'HETATM' lines in main part of file
* If providing cpp4 map files: All maps should be of the '.map' or '.cpp4' variety and labelled as:
    > ```
    > Aprot.pdb
    > Aprot_smiles.txt
    > Aprot_fofc.map
    > Aprot_2fofc.map
    > Aprot_event.cpp4
    > ```
* If you have multiple events, you can add use a `_event_0.cpp4` type naming convention and increment the number accordingly. Sorry! Here is an example of how you can do it though!
```bash
module load ccp4
cd /staging_dir/target/aligned
# The directory with the fragalysis stuff. ls should return folders for each ligand
folders=$(ls)
for i in $folders
  do
    str=$(echo ./$i/)
    maps=$(find $str -name *.map)
    ccp4s=$(find $str -name *.ccp4)
    for j in $maps
      do
        mapmask mapin $j mapout $j xyzin ./$i/$i.pdb << eof
        border 6
        end
        eof
      done

    for j in $ccp4s
      do
        mapmask mapin $j mapout $j xyzin ./$i/$i.pdb << eof
        border 6
        end
        eof
      done
  done
```

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
python fragalysis-api/fragalysis_api/xcimporter/xcimporter.py  -i [input directory] -o [output directory] -t [target name] -m 
```
Which will align the input files with respects to individual chains inside the .pdb files (`-m`) and save the output
to a folder specified by `-t [target name]`.

A description of the command line arguments are as follows:
- `-i`, `--in_dir` :  Input Directory
- `-o`, `--out_dir` : Output Directory
- `-v`, `--validate`: (Optional) Validate the Inputs
- `-m`, `--monomerize`: (Optional) Split the input PDBs into separate chains. E.G If a pdb has A and B chains it will create files pdb_name_A.pdb and pdb_name_B.pdb
- `-t`, `--target`: The name of the output folder to be saved in Output directory
- `-md`, `--metadata`: (Optional) Automatically populated a metadata.csv in the output directory to be fill in.
- `-b`, `--biomol`: (Optional) File path to plain text file that contains an optional header that you would like to be added to PDB files.
- `-r`, `--reference`: (Optional) The name/filepath of the pdb file you which to use as reference (can be PDB ID)
- `-c`, `--covalent`: (Optional) Handle Covalent attachments by extending output .mol file to include covalent attachment atoms. Requires modified smiles strings.

The terminal will let you know when the conversion has been successful and if there are any files that have been found to be incompatible with the API. We are working to minimize any incompatibilities.

The expected output of the xcimporter is a folder located at `[outputdirectory]/[target name]` which will
contain two main folders `aligned` and `crystallographic`. The `aligned` folder will contain multiple sub-folders, one per ligand. These will be labelled as the name of the `pdbfile_[lig_number]` or `pdbfile_[lig_number][chain_letter]` if monomer mode was specified.
Inside each of these subfolders contains all the information that is required by the fragalysis platform to display and elaborate a given ligand.
The `crystallographic` folder contains an exact copy of the data supplied to the above command.

#### 3.1 How to submit single PDB files for conversion to a fragalysis friendly format
If you would like to add individual pdb files to the results of the fragalysis api we can use:
```
python fragalysis-api/fragalysis_api/xcimporter/single_import.py --in_file=[pdbtobealigned.pdb] --out_dir=[output directory] --target [targetname] -m --reference[output directory/targetname/reference.pdb]
```
Where `in_file` is the filepath of a pdb file we would like to align to the rest of our results (with `_smiles.txt` and other files to be aligned in the same directory)
and `reference` is the pdb file that you would like to align files associated to `in_file`. If you had previous run xcimporter then a reference.pdb file should be located at `output direct/targetname/reference.pdb`
If `--reference` or `-r` are not specified then `output direct/targetname/reference.pdb` will be used according to what arguments you have specified.

You should be able to use `single_import` without having to delete any preexisting files. As the new outputs should overwrite what previously exists. Nice.

A description of the command line arguments for `single_import.py` are as follows: 
- `-i`, `--in_file` :  Input File
- `-o`, `--out_dir` : Output Directory
- `-v`, `--validate`: (Optional) Validate the Inputs
- `-m`, `--monomerize`: (Optional) Split the input PDBs into separate chains. E.G If a pdb has A and B chains it will create files pdb_name_A.pdb and pdb_name_B.pdb
- `-t`, `--target`: The name of the output folder to be saved in Output directory
- `-md`, `--metadata`: (Optional) Automatically populated a metadata.csv in the output directory to be fill in.
- `-b`, `--biomol`: (Optional) File path to plain text file that contains an optional header that you would like to be added to PDB files.
- `-r`, `--reference`: (Optional) The name/filepath of the pdb file you which to use as reference (can be PDB ID)
- `-c`, `--covalent`: (Optional) Handle Covalent attachments by extending output .mol file to include covalent attachment atoms. Requires modified smiles strings.
- `-sr`, `--self_reference`: (Optional) Indicate whether you want pdb files to align to themselves (for testing purposes)

#### 3.2 Running The Fragalysis API without alignment 
This is only available to the `single_import.py` method
If for whatever reason you decide that you would like to simply split a pdb file without the need of an alignment step you can use the `-sr` or `--selfreference` flags when using the fragalysis API. Or specify the `-r` or `--reference` to be itself. 

An example in bash to process a folder of pbds without aligned in bash would be:
```bash
$input=/path/to/input/folder
pdblist=$(ls input)
for pdb in $pdblist
do
  python fragalysis-api/fragalysis_api/xcimporter/single_import.py --in_file=$pdb --out_dir=[output directory] --target [targetname] -m --selfreference
done
```

Notably the `-sr` flag will take precedence over the `-r` flag.


### Preparing metadata for fragalysis upload
Each fragalysis upload requires a .csv file (called `metadata.csv`) to be stored inside the output `aligned`. It should have the following structure:


|<empty>|crystal_name|RealCrystalName|smiles|new_smiles|alternate_name|site_name|pdb_entry|
|---|---|---|---|---|---|---|---|
|1|prot-x0001_0A|prot-x0001|CCCC|<empty>|<empty>|Binding Site|<empty>|
|...|...|...|...|...|...|...|...|

This file will be generated by the api if using the `-md` flag or otherwise can be constructed manually.
- [First Column]: This doesn't have a column name, but just indicates a rownumber, fill 1 to n  
- __crystal_name__ (Required): The name of the crystal as it appear inside the aligned folder. e.g. prot-x0001_0A or prot-x0001_1A
- __RealCrystalName__ (Required): The name of the crystal that the ligand is derived from. e.g. prot-x0001, prot-x0002
- __smiles__ (Required): The SMILES string for the crystal as described in the `aligned/crystal_name/crystal_name_smiles.txt` file
- __new_smiles__ (can be left empty): If the SMILES string is slightly wrong or you want to specify a specific enantiomer to render, declare the new SMILES string here
- __alternate_name__ (can be left empty): If you have a particular title that you would like to render alongside the crystal_name
- __site_name__ (Recommended): This is the most important field. All rows in this column need a value! e.g. Binding Site. So that Fragalysis will know where to cluster the fragments together. This will be rendered alphabetically so if you want to preserved a specific order of sites, prefix the site labels with alphanumeric combinations e.g. "A1 - Binding Site", "B2 - Allosteric Site". If entire column is left empty then fragalysis will cluster ligands together by center of mass and will use uninformative site labels.
- __pdb_entry__ (can be left empty): The PDB code relating to the Crystal.

The rendered table line in plain text will look like this:
```
,crystal_name,RealCrystalName,smiles,new_smiles,alternate_name,site_name,pdb_entry
1,prot-x0001_0A,prot-x0001,CCCC,,,Binding Site,
```
### Who we are
We are the Fragment 5, a group of students at the University of Oxford.

Anna :whale:   
Maranga :fire:  
George a.k.a Joff Boff :poop:  
Tobias :cow:  
Alister :panda_face:  

We are looked after by Rachael :crown:. 

Finally, the code quality was ruined by Tyler.




