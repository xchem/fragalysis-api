import os
import argparse
from fragalysis_preproc.data import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                     "Get apo pdbs and molecule sdfs in csv format for an input target"
                                     " from https://fragalyysis.diamond.ac.uk")
    parser.add_argument('-t', '--target', help='Name of target on Fragalysis')
    parser.add_argument('-o', '--output', help='Output csv file')

    parser.add_argument('-s', '--save', action='store_true')
    parser.add_argument('-d', '--directory', help='directory to save pdbs and sdfs')

    args = parser.parse_args()

    if args.target:
        target = args.target
    else:
        raise Exception('No target input provided - please see help! (--help)')

    if not args.output:
        raise Exception('No output csv provided - please see help! (--help)')

    search = GetMoleculesData()
    search.get_target_ids(target)
    search.get_all_mol_responses()
    dct = search.convert_mols_to_dict()

    if args.save:
        print("\nResults will be saved into the directory passed by -d (--directory)\n"
              "with a sub-folder for each hit, containing:\n"
              "1. apo pdb file - <fragalysis_code>.pdb\n"
              "2. sdf of bound ligand - <fragalysis_code>.sdf\n\n"
              "A csv file containing the smiles string for each code will be saved\n"
              "in the directory passed by -d (--directory)\n")
        if not args.directory:
            raise Exception('please specify an output directory!')
        if args.directory:
            if not os.path.isdir(args.directory):
                os.mkdir(args.directory)
            for row in dct.itertuples():
                path = os.path.join(args.directory, row.code)
                if not os.path.isdir(path):
                    os.mkdir(path)
                pdb_path = os.path.join(path, str(row.code + '.pdb'))
                with open(pdb_path, 'w') as f:
                    f.write(row.pdb)
                sdf_path = os.path.join(path, str(row.code + '.sdf'))
                with open(sdf_path, 'w') as f:
                    f.write(row.sdf)
            dct.to_csv(os.path.join(args.directory, args.output.split('/')[-1]), columns=['code', 'smiles'])
    if not args.save:
        dct.to_csv(args.output)


