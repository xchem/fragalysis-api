import argparse
from sys import exit
import os

from shutil import copyfile
from fragalysis_api import Validate, Align2, Monomerize
from fragalysis_api import set_up

from distutils.dir_util import copy_tree


def xcimporter(in_dir, out_dir, target, metadata=False, validate=False, monomerize=False):
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

    if validate:
        validation = Validate(in_dir)

        if not bool(validation.is_pdbs_valid):
            print("Input files are invalid!!")
            exit

        if not validation.does_dir_exist:
            print("Input dir doesn't exist.")
            exit

        if not validation.is_there_a_pdb_in_dir:
            print("No PDB file in input dir.")
            exit

    print("Making output directories")
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
        os.makedirs(os.path.join(out_dir, "tmp"))

    in_dir2 = in_dir

    if monomerize:
        print("Monomerizing input structures")
        out = os.path.join(out_dir, 'mono/')
        if not os.path.isdir(out):
            os.makedirs(out)
        mono = Monomerize(directory=in_dir, outdir=out)
        mono.monomerize_all()
        in_dir = out

    pdb_smiles_dict = {'pdb': [], 'smiles': []}

    for f in os.listdir(in_dir):
        if '.pdb' in f:
            pdb_smiles_dict['pdb'].append(os.path.join(in_dir, f))
            print(os.path.join(in_dir, f).replace('.pdb', '_smiles.txt'))
            if os.path.isfile(os.path.join(in_dir, f).replace('.pdb', '_smiles.txt')):
                pdb_smiles_dict['smiles'].append(os.path.join(in_dir, f).replace('.pdb', '_smiles.txt'))
            else:
                pdb_smiles_dict['smiles'].append(None)

    print(pdb_smiles_dict['smiles'])
    print("Aligning protein structures")
    print('New Stuff')
    structure = Align2(in_dir, "", monomerize)
    structure.align(out_dir=os.path.join(out_dir, "tmp"))

    for smiles_file in pdb_smiles_dict['smiles']:
        if smiles_file:
            print(smiles_file)
            copyfile(smiles_file, os.path.join(os.path.join(out_dir, "tmp", smiles_file.split('/')[-1])))
            print(os.path.join(out_dir, "tmp", smiles_file.split('/')[-1]))

    aligned_dict = {'bound_pdb':[], 'smiles':[]}

    for f in os.listdir(os.path.join(out_dir, "tmp")):
        if '.pdb' in f:
            aligned_dict['bound_pdb'].append(os.path.join(out_dir, "tmp",f))
            if os.path.isfile(os.path.join(out_dir, "tmp",f).replace('_bound.pdb', '_smiles.txt')):
                aligned_dict['smiles'].append(os.path.join(out_dir, "tmp",f).replace('_bound.pdb', '_smiles.txt'))
            else:
                aligned_dict['smiles'].append(None)


    print(aligned_dict['smiles'])

    # Moving Cut will now be moved to outside of api
    #print("Cutting Maps around ligands")
    #out = os.path.join(out_dir, 'cut/')
    #if not os.path.isdir(out):
    #    os.makedirs(out)
    #cutmaps = CutMaps(in_dir=in_dir, out_dir=out, monomerize=monomerize)
    #cutmaps.cut_maps()

    print("Identifying ligands")
    for aligned, smiles in list(zip(aligned_dict['bound_pdb'], aligned_dict['smiles'])):
        try:
            if smiles:
                _ = set_up(target_name=target, infile=os.path.abspath(aligned),
                             out_dir=out_dir, monomerize=monomerize, smiles_file=os.path.abspath(smiles))
                
            else:
                _ = set_up(target_name=target, infile=os.path.abspath(aligned),
                out_dir=out_dir, monomerize=monomerize)
                
        except AssertionError:
            print(aligned, "is not suitable, please consider removal or editing")
            for file in os.listdir(os.path.join(out_dir, "tmp")):
                if str(aligned) in file:
                    os.remove(os.path.join(out_dir, "tmp", str(file)))

    if metadata == '1':
        print("Preparing metadata file")
        metadata_fp = os.path.join(out_dir, target, "metadata.csv")
        path = os.path.join(out_dir, target)

        with open(metadata_fp, 'w+') as f:
            for root, dirs_list, files_list in os.walk(path):
                for file_name in files_list:
                    if file_name.endswith('meta.csv'):
                        csv_path = os.path.join(root, file_name)
                        for line in open(csv_path, 'r'):
                            f.write(line)

    # Move input files into Target/crystallographic folder
    copy_tree(in_dir2, os.path.join(out_dir, target, 'crystallographic') )

    print("Files are now in a fragalysis friendly format!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--in_dir",
        default=os.path.join("..", "..", "data", "xcimporter", "input"),
        help="Input directory",
        required=True,
    )
    parser.add_argument(
        "-o",
        "--out_dir",
        default=os.path.join("..", "..", "data", "xcimporter", "output"),
        help="Output directory",
        required=True,
    )
    parser.add_argument(
        "-v", "--validate", action="store_true", default=False, help="Validate input"
    )
    parser.add_argument(
        "-m", "--monomerize", action="store_true", default=False, help="Monomerize input"
    )
    parser.add_argument("-t", "--target", help="Target name", required=True)
    parser.add_argument("-md", "--metadata", help="Metadata output", default=False)

    args = vars(parser.parse_args())

    # user_id = args['user_id']
    in_dir = args["in_dir"]
    out_dir = args["out_dir"]
    validate = args["validate"]
    monomerize = args["monomerize"]
    target = args["target"]
    metadata = args["metadata"]

    if in_dir == os.path.join("..", "..", "data", "xcimporter", "input"):
        print("Using the default input directory ", in_dir)
    if out_dir == os.path.join("..", "..", "data", "xcimporter", "output"):
        print("Using the default input directory ", out_dir)

    xcimporter(in_dir=in_dir, out_dir=out_dir, target=target, validate=validate, monomerize=monomerize, metadata=metadata)

    fix_pdb = open(os.path.join(out_dir, target, 'aligned', 'pdb_file_failures.txt'), 'w')

    for target_file in os.listdir(os.path.join(out_dir, target)):
        if target_file != 'pdb_file_failures.txt' and len(os.listdir(os.path.join(out_dir, target, target_file))) < 2:
            rmtree(os.path.join(out_dir, target, target_file))
            fix_pdb.write(target_file.split('-')[1]+'\n')

    fix_pdb.close()
    print('For files that we were unable to process, look at the pdb_file_failures.txt file in your results directory.'
          ' These files were unable to produce RDKit molecules, so the error likely lies in the way the ligand atoms or'
          'the conect files have been written in the pdb file')
