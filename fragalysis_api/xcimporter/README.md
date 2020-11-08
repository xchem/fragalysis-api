Preparing files for fragalysis
------------------------------
The script ``fragalysis-api/fragalysis_api/xcimporter/xcimporter.py`` can be used to prepare bound-state pdb files for upload
to fragalysis. 

Usage:
------
options:  

```python:

-i, --in_dir: input directory where bound-state pdb files are
-o, --out_dir: output directory where target folder containing data should be written to
-m, --monomerize: whether to monomorize all input structures to align everything to the same reference frame
-t, --target: the name of the target, as it is to appear on the fragalysis landing page
-md, --metadata: whether to prepare a metadata template file (1=True, 0=False)
```

The input directory should contain a bound state pdb file for eacch crystal structure of the target, and a ``.txt`` file of 
the same filename for each structure containing the smiles string of the bound ligand.  

For example:

```
Mpro
- Mpro-x123.pdb
- Mpro-x123.txt
```

Example to run:

```
python fragalysis-api/fragalysis_api/xcimporter/xcimporter.py --in_dir=/dls/science/groups/i04-1/fragprep/input/Mpro/ --out_dir=/dls/science/groups/i04-1/fragprep/staging --target Mpro
```

After running, the data will be output in a directory named by the ``target`` option, under the location specified by ``out_dir``
