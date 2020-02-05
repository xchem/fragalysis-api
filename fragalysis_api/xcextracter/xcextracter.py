from getdata import GetTargetsData, GetPdbData, GetMoleculesData
from graph import GraphRequest


def xcextracter(target_name):
    import pandas as pd

    search = GetMoleculesData()

    # get the target id's for further operations for NUDT5A
    search.set_target_id(target=target_name)

    # set the molecule url for the current instance of GetMoleculesData to the first id in the id_list from above
    search.set_molecule_url()

    # get a json response from the url we set above
    search.set_mol_data()

    # get all molecule data (json) related to all of the ids in id_list
    search.set_complete_mol_data()

    a_df = pd.DataFrame(search.get_complete_mol_data)
    print(a_df)

    return a_df


def xcgraphcreator(target_smiles):

    search = GraphRequest()

    # set the search smiles
    search.set_smiles_url(smiles=target_smiles)

    results = search.get_graph_json()

    return results


if '__main__' == __name__:

    #xcextracter(target_name='ATAD')

    xcgraphcreator(target_smiles='O=C(Nc1ccccc1)Nc1cccnc1')