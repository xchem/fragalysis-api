# Pulling Data From Fragalysis

## Settings
Saved in config.ini in top level - default in this repo (go up a folder from here!)

**[fragalysis]** - info about fragalysis deployment to search against

url = url for deployment to search against. Default = fragalysis.diamond.ac.uk

**[graph]** - parameters for searching the graph network

search = restful endpoint (extension of [fragalysis] url) to query the graph network

query = search string in url format for search field to query against. Currently a smiles string

**[targets]** - parameters for querying target (protein) information

search = restful endpoint (extension of [fragalysis] url) to query the fragalysis db for a list of targets

query = search string in url format for search field to query against. Currently the target name

**[molecules]** - parameters for querying molecule (ligand) information

search = search = restful endpoint (extension of [fragalysis] url) to query the fragalysis db for molecule info

query = search string in url format for search field to query against. Currently the target id (to get all related molecules)

**[pdb]** - parameters for getting pdb files

search = search = restful endpoint (extension of [fragalysis] url) to query the fragalysis db for pdb files

query = search string in url format for search field to query against. Currently the code of the target molecule - usually similar to the crystal name from an experiment

## Structural data - apo structure and sdf for each hit in a named target 
### Using the API

**1. Getting target data**

First, init an instance of GetTargetsData from fragalysis_preproc/data.py:

```
from data import GetTargetsData

search = GetTargetsData()
```

Next, we need to set a target name that we want to search for. This is done by setting the target_name_url attribute:

```
# e.g. set NUDT5A as the target for the search
search.set_target_name_url(target='NUDT5A')
```

Now, we can get a json response, which gives back the info specified by the fragalysis RESTful api (see <insert documentation link when it exists>):

```
# get the json response
results = search.get_target_json()
```

Finally, we can parse the json with a built-in function to return id's (db id) for the entries relating to the target we searched. We can use this to make further queries to find data in other tables relating to this target:

```
# get a list of integers corresponding to the pk's for all entries related to the search target
id_list = search.get_target_id_list()
```

**2. Getting pdb data - currently apo pdb for given crystal name**

First, init an instance of GetPdbData from fragalysis_preproc/data.py

```
from data import GetPdbData

search = GetPdbData()
```

Next, we can get a string representation of the PDB file corresponding to the apo structure for the bound-state structure, given an input 'code', i.e. crystal name

```
# get apo pdb for bound-state of NUDT5A-x0123 in string representation
apo_pdb = search.get_pdb_file(code='NUDT5A-x0123')
```

**3. Getting molecule data - all data relating to molecule from molecule table**

First, init an instance of GetMoleculesData from fragalysis_preproc/data.py

```
from data import GetMoleculesData

search = GetMoleculesData()
```

Further searches require a target id, which can be found with the api, as shown in point 1. Alternativley, this can be done from the instance we just set up, for a given target name:

```
# get the target id's for further operations for NUDT5A
id_list = search.get_target_ids(target='NUDT5A')
```

The GetMoleculesData object can either operate on one object at a time, or the entire id list.

Here, I'll show an example that will just cover the first target id in `id_list`

To search the API for moleaules relating to the target id, we need to sett the molecule url:

```
# set the molecule url for the current instance of GetMoleculesData to the first id in the id_list from above
url = search.set_molecule_url(target_id=id_list[0])
```

Now, we can get the data held in the Molecules table in the fragalysis db in json format (we need to use the returned url from the above snippet. The reason why will become clear later)

```
# get a json response from the url we set above
results = get_molecules_json(url=url)
```

If we want to operate on all target id's related to a target, we can use a built-in function. Make sure to first set the target_ids parameter with `get_target_ids` as shown earlier.

```
# get all molecule data (json) related to all of the ids in id_list
results_table = search.get_all_mol_responses()
```

This returns a pandas dataframe containing: code(crystal name), pdb(apo pdb file as string), sdf(sdf of the bound ligand as string) and smiles(smiles string of the bound ligand)

You can then work on this dataframe however you wish. 

Alternativley, there is a script which can be used to download the apo pdb's and sdf files, and a csv file containing the smiles strings for each hit for a given target, as shown in the next section

### Using script to download

The script lives in the top level, and is called get_target_data.py

You can see how to run this script with `python get_target_data.py --help`

## Graph Network search <in progress>
We can search the graph network currently deployed at the fragalysis deployment specified in config.ini by smiles string

First, we initialise an instance of GraphRequest from fragalysis_preproc/graph.py

```
from graph import GraphRequest

search = GraphRequest()
```

Next, we set the smiles string we want to search

```
# set the search smiles
search.set_smiles_url(smiles='O=C(Nc1ccccc1)Nc1cccnc1')
```

Next, we can get the response from the graph network in json format

```
results = search.get_graph_json()

TODO: convert json reponse from graph into a list of smiles strings for follow up compounds


