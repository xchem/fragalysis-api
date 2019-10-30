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