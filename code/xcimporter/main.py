from align import Align
from conversion_pdb_mol import set_up
from validate import Validate
import os

if __name__ == "__main__":
    user_id = input("User ID:")
    in_dir = '../../data/xcimporter/input/'
    out_dir = '../../data/xcimporter/output/'
    pdb_list = ['6epu', '6epv', '6epx', '6hi3']

    if not os.path.exists(os.path.join(out_dir, str(user_id), '/')):
        os.makedirs(os.path.join(out_dir, str(user_id), '/'))
        os.makedirs(os.path.join(out_dir, str(user_id), 'tmp/'))

    validation = Validate(in_dir)
    validation.validate_pdbs

    struc = Align(os.path.join(in_dir, str(user_id)), pdb_ref='')
    struc.align(os.path.join(out_dir, str(user_id), 'tmp/'))

    for i in pdb_list:
        new = set_up(i, str(user_id))
