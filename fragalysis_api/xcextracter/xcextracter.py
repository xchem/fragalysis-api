from fragalysis_api import GetMoleculesData, can_connect
import pandas as pd
import sys


def xcextracter(target_name):

    #checking the fragalysis website is live
    if can_connect():
        print('Can connect to Fragalysis')
    else:
        print('Cannot connect to Fragalysis')
        sys.exit()

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

    return a_df


if '__main__' == __name__:

    results = xcextracter(target_name='NUDT5A')
    print(results)

