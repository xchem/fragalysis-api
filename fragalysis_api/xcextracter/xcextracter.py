from getdata import GetTargetsData, GetPdbData, GetMoleculesData
from graph import GraphRequest


def xcextracter(target_name):

    # Initiate the object
    search = GetTargetsData()
    # Set the target into the object
    search.set_target_name_url(target=target_name)
    # get the json response (it internally calls the frag database)
    results = search.get_target_json()

    print(results)

    # get a list of integers corresponding to the pk's for all entries related to the search target
    id_list = search.get_target_id_list()

    print(id_list)

    search = GetPdbData()

    # get apo pdb for bound-state of NUDT5A-x0123 in string representation
    apo_pdb = search.get_pdb_file(code='NUDT7A-x0140_1')

    search = GetMoleculesData()

    # get the target id's for further operations for NUDT5A
    id_list = search.get_target_ids(target='NUDT5A')

    # set the molecule url for the current instance of GetMoleculesData to the first id in the id_list from above
    url = search.set_molecule_url(target_id=id_list[0])

    # get a json response from the url we set above
    results = search.get_molecules_json(url=url)

    # get all molecule data (json) related to all of the ids in id_list
    results_table = search.get_all_mol_responses()

    print(results_table)

    search = GraphRequest()

    # set the search smiles
    search.set_smiles_url(smiles='O=C(Nc1ccccc1)Nc1cccnc1')

    results = search.get_graph_json()

    import pandas as pd

    return results



if '__main__' == __name__:

    xcextracter(target_name='ATAD')