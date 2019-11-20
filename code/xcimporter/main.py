from align import Align
from conversion_pdb_mol import set_up
from validate import Validate
from create_directory import create_directory
import os
from shutil import rmtree

if __name__ == "__main__":
    user_id = input("User ID:")
    in_dir = '../../data/xcimporter/input/'
    out_dir = '../../data/xcimporter/output/'
    pdb_list = []

    for file in os.listdir(os.path.join(in_dir, str(user_id))):
        pdb_list.append(file[:-4])

    if not os.path.exists(os.path.join(out_dir, str(user_id), '/')):
        os.makedirs(os.path.join(out_dir, str(user_id), '/'))
        os.makedirs(os.path.join(out_dir, str(user_id), 'tmp/'))

    validation = Validate(os.path.join(in_dir, str(user_id)))
    validation.validate_pdbs

    struc = Align(os.path.join(in_dir, str(user_id)), pdb_ref='')
    struc.align(os.path.join(out_dir, str(user_id), 'tmp/'))

    for i in pdb_list:
        new = set_up(i, str(user_id))

    create_directory(str(user_id), os.path.join(out_dir, str(user_id), 'tmp/'))

    rmtree(os.path.join(out_dir, str(user_id), 'tmp/'))
