import os
import argparse
from Bio.PDB import PDBList


def pdb_importer(data_dir, user_id, pdb_code):
    if not os.path.exists(os.path.join(data_dir, user_id)):
        os.mkdir(os.path.join(data_dir, user_id))
        print('making directory')

    if not os.path.exists(os.path.join(data_dir, user_id, pdb_code + '.pdb')):
        pdbl = PDBList()
        pdbl.retrieve_pdb_file(pdb_code, pdir=os.path.join(data_dir, user_id), file_format='pdb')
        os.rename(os.path.join(data_dir, user_id, 'pdb' + pdb_code + '.ent'),
                  os.path.join(data_dir, user_id, pdb_code + '.pdb'))
    else:
        print('File is already downloaded')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-id', '--user_id', required=True, help='user ID')
    parser.add_argument('-pdb', '--pdb', required=True, help='pdb code you want to import')
    args = vars(parser.parse_args())

    user_id = args['user_id']
    data_dir = os.path.join('..', '..', 'data', 'xcimporter', 'input')
    pdb_code = args['pdb']

    pdb_importer(data_dir, user_id, pdb_code)
