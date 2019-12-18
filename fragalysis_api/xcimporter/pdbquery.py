from pypdb import get_blast2, get_ligands
import argparse
import json
import os
from pdbimporter import pdb_importer


def get_similar(pdb_code, chain_id):
    """
    queries the pdb for other files with similar protein structure but different ligands bound

    params: pdb code, chain id
    returns: prints list of pdb codes and ligands of similar protein structures, or saves this as a dictionary in a json file
    """
    similar = [pdb_code]
    result = get_blast2(pdb_code, chain_id=chain_id)
    codes = result[0]
    info = result[1]
    for i in range(len(codes)):
        if '100%' in info[i].text:
            similar.append(codes[i])

    similar_with_ligs = {}

    for i in similar:
        ligands = get_ligands(i)  # finds ligands bound to the protein in that particular pdb file
        try:
            if ligands['ligandInfo']['ligand'][0]['@chemicalID'] not in similar_with_ligs:  # if ligand not already seen
                similar_with_ligs[ligands['ligandInfo']['ligand'][0]['@chemicalID']] = [i.lower()]   # save pdb code
            else:
                similar_with_ligs[ligands['ligandInfo']['ligand'][0]['@chemicalID']].append(i.lower())
        except KeyError:
            if ligands['ligandInfo']['ligand']['@chemicalID'] not in similar_with_ligs:
                similar_with_ligs[ligands['ligandInfo']['ligand']['@chemicalID']] = [i.lower()]
            else:
                similar_with_ligs[ligands['ligandInfo']['ligand']['@chemicalID']].append(i.lower())
        except TypeError:
            pass

    for i in similar_with_ligs:
        similar_with_ligs[i] = list(set(similar_with_ligs[i]))
    print(len(similar_with_ligs)-1, 'different ligands have been found to bind to this protein')
    print_list = input('would you like to view these ligands and the pdb codes? y/n : ')
    if print_list == 'y' or print_list == 'Y' or print_list == 'yes':
        for i in similar_with_ligs:
            print(i, similar_with_ligs[i])
    save_file = input('would you like to save this list in a json file? y/n : ')
    if save_file == 'y' or save_file == 'Y' or save_file == 'yes':
        user_id = input('user id: ')
        if not os.path.exists(os.path.join('..', '..', 'data', 'xcimporter', 'other', user_id)):
            os.mkdir(os.path.join('..', '..', 'data', 'xcimporter', 'other', user_id))
        json.dump(similar_with_ligs, open(os.path.join('..', '..', 'data', 'xcimporter', 'other', user_id, pdb_code+'.json'), 'w'))
    import_pdb = input('would you like to import these pdb files? y/n : ')
    if import_pdb == 'y' or import_pdb == 'Y' or import_pdb == 'yes':
        user_id = input('user id: ')
        data_dir = os.path.join('..', '..', 'data', 'xcimporter', 'input')
        for i in similar_with_ligs:
            for j in similar_with_ligs[i]:
                pdb_importer(data_dir, user_id, j)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb', '--pdb', required=True, help='pdb code to search for similar')
    parser.add_argument('-chain', '--chain', required=True, help='chain of structure you want to look at')
    args = vars(parser.parse_args())

    pdb_code = args['pdb']
    chain_id = args['chain']

    get_similar(pdb_code, chain_id)
