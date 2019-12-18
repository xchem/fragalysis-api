from pypdb import get_blast2, get_ligands
import argparse
import json
import os
from pdbimporter import ImportPdb

class Query:
    def __init__(self, pdb, chain):
        self.pdb_code = pdb
        self.chain_id = chain
        self.similar = [self.pdb_code]
        self.match_ligs = {}


    def get_matching_proteins(self):
        """
        queries the pdb for other files with similar protein structure but different ligands bound

        params: pdb code, chain id
        returns: prints list of pdb codes and ligands of similar protein structures, or saves this as a dictionary in a json file
        """
        result = get_blast2(self.pdb_code, chain_id=self.chain_id)
        codes = result[0]
        info = result[1]
        for i in range(len(codes)):
            if '100%' in info[i].text:
                self.similar.append(codes[i])

    def get_ligands(self):

        for i in self.similar:
            ligands = get_ligands(i)  # finds ligands bound to the protein in that particular pdb file
            try:
                if ligands['ligandInfo']['ligand'][0]['@chemicalID'] not in self.match_ligs:  # if ligand not already seen
                    self.match_ligs[ligands['ligandInfo']['ligand'][0]['@chemicalID']] = [i.lower()]   # save pdb code
                else:
                    self.match_ligs[ligands['ligandInfo']['ligand'][0]['@chemicalID']].append(i.lower())
            except KeyError:
                if ligands['ligandInfo']['ligand']['@chemicalID'] not in self.match_ligs:
                    self.match_ligs[ligands['ligandInfo']['ligand']['@chemicalID']] = [i.lower()]
                else:
                    self.match_ligs[ligands['ligandInfo']['ligand']['@chemicalID']].append(i.lower())
            except TypeError:
                pass

    def print_number_ligs(self):
        for i in self.match_ligs:
            self.match_ligs[i] = list(set(self.match_ligs[i]))
        print(len(self.match_ligs)-1, 'different ligands have been found to bind to this protein')
        return len(self.match_ligs)

    def view_ligands(self):
        for i in self.match_ligs:
            print(i, self.match_ligs[i])
        return len(self.match_ligs)

    def save_dictionary(self, user):
        if not os.path.exists(os.path.join('..', '..', 'data', 'xcimporter', 'other', user)):
            os.mkdir(os.path.join('..', '..', 'data', 'xcimporter', 'other', user))
        json.dump(self.match_ligs, open(os.path.join('..', '..', 'data', 'xcimporter', 'other', user_id, self.pdb_code+'.json'), 'w'))

    def import_pdbs(self, user):
        data_dir = os.path.join('..', '..', 'data', 'xcimporter', 'input')
        for i in self.match_ligs:
            for j in self.match_ligs[i]:
                initial = ImportPdb(data_dir, user, j)
                initial.pdb_importer()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-pdb', '--pdb', required=True, help='pdb code to search for similar')
    parser.add_argument('-chain', '--chain', required=True, help='chain of structure you want to look at')
    args = vars(parser.parse_args())

    pdb_code = args['pdb']
    chain_id = args['chain']

    query_obj = Query(pdb_code, chain_id)
    query_obj.get_matching_proteins()
    query_obj.get_ligands()
    query_obj.print_number_ligs()
    print_list = input('would you like to view these ligands and the pdb codes? y/n : ')
    if print_list == 'y' or print_list == 'Y' or print_list == 'yes':
        query_obj.view_ligands()
    save_file = input('would you like to save this list in a json file? y/n : ')
    if save_file == 'y' or save_file == 'Y' or save_file == 'yes':
        user_id = input('user id: ')
        query_obj.save_dictionary(user_id)
    import_pdb = input('would you like to import these pdb files? y/n : ')
    if import_pdb == 'y' or import_pdb == 'Y' or import_pdb == 'yes':
        user_id = input('user id: ')
        query_obj.import_pdbs(user_id)
