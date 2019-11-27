from align import Align
from conversion_pdb_mol import set_up
from validate import Validate
from create_directory import create_directory
import os
from shutil import rmtree
import argparse

if __name__ == "__main__":
    #make these commandline options using argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-id", "--user_id", help="specify user ID", required = True)
    parser.add_argument("-i", "--in_dir", default='../../data/xcimporter/input/', help="specify the directory the PDB files is in", required = False)
    parser.add_argument("-o", "--out_dir", default='../../data/xcimporter/output/', help="specify the directory to output the xchem format PDBs to", required = False)
    args = parser.parse_args()
    if args.in_dir == '../../data/xcimporter/input/':
        print('Using the default input directory ', args.in_dir)
    if args.out_dir == '../../data/xcimporter/output/':
        print('Using the default input directory ', args.out_dir)
    pdb_list = []

    for file in os.listdir(os.path.join(args.in_dir, str(args.user_id))):
        pdb_list.append(file[:-4])

    if not os.path.exists(os.path.join(args.out_dir, str(args.user_id), '/')):
        os.makedirs(os.path.join(args.out_dir, str(args.user_id), '/'))
        os.makedirs(os.path.join(args.out_dir, str(args.user_id), 'tmp/'))

    validation = Validate(os.path.join(args.in_dir, str(args.user_id)))
    validation.validate_pdbs

    struc = Align(os.path.join(args.in_dir, str(args.user_id)), pdb_ref='')
    struc.align(os.path.join(args.out_dir, str(args.user_id), 'tmp/'))

    for i in pdb_list:
        new = set_up(i, str(args.user_id))

    create_directory(str(args.user_id), os.path.join(args.out_dir, str(args.user_id), 'tmp/'))

    rmtree(os.path.join(args.out_dir, str(args.user_id), 'tmp/'))
