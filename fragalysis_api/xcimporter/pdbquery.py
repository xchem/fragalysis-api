from pypdb import get_blast2, get_ligands
import argparse
import json
import os


def get_similar(pdb_code, chain_id):
    similar_with_ligs = []
    ligs = []
    similar = [pdb_code]
    similar += get_blast2(pdb_code, chain_id=chain_id)[0]

    for i in similar[:20]:
        ligands = get_ligands(i)
        try:
            if ligands['ligandInfo']['ligand'][0]['@chemicalID'] not in ligs:
                similar_with_ligs.append(i)
                ligs.append(ligands['ligandInfo']['ligand'][0]['@chemicalID'])
        except KeyError:
            if ligands['ligandInfo']['ligand']['@chemicalID'] not in ligs:
                similar_with_ligs.append(i)
                ligs.append(ligands['ligandInfo']['ligand']['@chemicalID'])
        except TypeError:
            pass

    print(len(similar_with_ligs)-1, 'pdb files with the same structure and a different ligand bound have been found')
    print_list = input('would you like to view these pdb codes and their ligands? y/n : ')
    if print_list == 'y' or print_list == 'Y' or print_list == 'yes':
        for i in range(len(similar_with_ligs)):
            print(similar_with_ligs[i], ligs[i])
    save_file = input('would you like to save this list in a json file? y/n : ')
    if save_file == 'y' or save_file == 'Y' or save_file == 'yes':
        user_id = input('user id: ')
        if not os.path.exists(os.path.join('..', '..', 'data', 'xcimporter', 'other', user_id)):
            os.mkdir(os.path.join('..', '..', 'data', 'xcimporter', 'other', user_id))
        comb_list = {}
        for i in range(len(similar_with_ligs)):
            comb_list[similar_with_ligs[i]] = ligs[i]
        json.dump(comb_list, open(os.path.join('..', '..', 'data', 'xcimporter', 'other', user_id, pdb_code+'.json'), 'w'))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb', '--pdb', required=True, help='pdb code to search for similar')
    parser.add_argument('-chain', '--chain', required=True, help='chain of structure you want to look at')
    args = vars(parser.parse_args())

    pdb_code = args['pdb']
    chain_id = args['chain']

    get_similar(pdb_code, chain_id)
