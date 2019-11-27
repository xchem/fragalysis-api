from align import Align
from conversion_pdb_mol import set_up
from validate import Validate
from create_directory import create_directory
import os
from shutil import rmtree
import argparse


def xcimporter(user_id, in_dir, out_dir):
    """Formats a lists of PDB files into fragalysis friendly format.
    1. Validates the naming of the pdbs.
    2. It aligns the pdbs (_bound.pdb file).
    3. It cleans the pdbs (removes solvents and ions).
    4. Splits the pdbs into ligands (.mol, .sdf and .pdb formats) and protein files (_apo.pdb file).
    5. Orders all the files into the correct directory structure required for fragalysis.

    :param user_id: data ID given to the input. Should be present in in_dir.
    :param in_dir: Directory containing data ID directories.
    :param out_dir: Directory containing processed pdbs (will be created if it doesn't exists).
    :return:
    """
    validation = Validate(os.path.join(in_dir, str(user_id)))

    if bool(validation.is_pdbs_val):
        exit()


    pdb_list = []

    for file in os.listdir(os.path.join(in_dir, str(user_id))):
        pdb_list.append(file[:-4])

    if not os.path.exists(os.path.join(out_dir, str(user_id), '/')):
        os.makedirs(os.path.join(out_dir, str(user_id), '/'))
        os.makedirs(os.path.join(out_dir, str(user_id), 'tmp/'))

    struc = Align(os.path.join(in_dir, str(user_id)), pdb_ref='')
    struc.align(os.path.join(out_dir, str(user_id), 'tmp/'))

    for i in pdb_list:
        try:
            new = set_up(i, str(user_id))
        except AssertionError:
            print(i, 'is not suitable, please consider removal or editing')
            for file in os.listdir(os.path.join(out_dir, str(user_id), 'tmp/')):
                if str(i) in file:
                    os.remove(os.path.join(out_dir, str(user_id), 'tmp', str(file)))
            pass

    create_directory(str(user_id), os.path.join(out_dir, str(user_id), 'tmp/'))

    rmtree(os.path.join(out_dir, str(user_id), 'tmp/'))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-id', '--user_id', required=True, help='Description for foo argument')
    parser.add_argument('-i', '--in_dir', default='../../data/xcimporter/input/', help='Description for bar argument')
    parser.add_argument('-o', '--out_dir', default='../../data/xcimporter/output/', help='Description for bar argument')
    args = vars(parser.parse_args())
    
    user_id = args['user_id']
    in_dir = args['in_dir']
    out_dir = args['out_dir']

    if in_dir == '../../data/xcimporter/input/':
        print('Using the default input directory ', in_dir)
    if out_dir == '../../data/xcimporter/output/':
        print('Using the default input directory ', out_dir)

    xcimporter(user_id, in_dir, out_dir)