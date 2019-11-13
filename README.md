# fragalysis-api



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
#### :warning:**WARNING**:warning: *Currently aglin is only expecting 1 chain!!!!!!!!*:loudspeaker:

##:scroll:Enforced rules
* The pdb file shall not be greater than 5mb.
* The pdb filename shall not contain non English lanugague ascii characters and shall be between 4 and 20 characters in length.
* Each pdb file for alignment shall contain the same number of chains.
* All pdb files to be aligned must be in the same directory.
* If manually selecting a file for reference it must be in the same director as all pdb files for alignment. 


