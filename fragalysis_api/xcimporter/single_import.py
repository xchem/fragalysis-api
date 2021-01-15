import argparse
import glob
import os
import shutil

from fragalysis_api import Align, Monomerize, set_up


def import_single_file(in_file, out_dir, target, monomerize, reference, biomol=None, covalent=False, self_ref=False):

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    if not os.path.isdir(os.path.join(out_dir, f"tmp{target}/")):
        os.makedirs(os.path.join(out_dir, f"tmp{target}/"))

    in_dir = os.path.dirname(in_file)

    pdb_smiles_dict = {'pdb': [], 'smiles': []}

    if monomerize:
        print("Monomerizing input structure")
        out = os.path.join(out_dir, f'mono{target}/')
        if not os.path.isdir(out):
            os.makedirs(out)

        mono = Monomerize(directory=in_dir, outdir=out)
        mono.monomerize_single(file=in_file)
        in_dir = out
        for f in os.listdir(in_dir):
            if '.pdb' in f:
                pdb_smiles_dict['pdb'].append(os.path.join(in_dir, f))
                if os.path.isfile(os.path.join(in_dir, f).replace('.pdb', '_smiles.txt')):
                    pdb_smiles_dict['smiles'].append(os.path.join(in_dir, f).replace('.pdb', '_smiles.txt'))
                else:
                    pdb_smiles_dict['smiles'].append(None)
    else:
        pdb_smiles_dict['pdb'].append(in_file)
        if os.path.isfile(in_file.replace('.pdb', '_smiles.txt')):
            pdb_smiles_dict['smiles'].append(in_file.replace('.pdb', '_smiles.txt'))
        else:
            pdb_smiles_dict['smiles'].append(None)

    print("Aligning to Reference")
    structure = Align(in_dir, "", monomerize)

    for i in pdb_smiles_dict['pdb']:
        if self_ref:
            structure.align_to_reference(i, i, out_dir=os.path.join(out_dir, f"tmp{target}"))
        else:
            structure.align_to_reference(i, reference, out_dir=os.path.join(out_dir, f"tmp{target}"))

    for smiles_file in pdb_smiles_dict['smiles']:
        if smiles_file:
            print(smiles_file)
            shutil.copyfile(smiles_file, os.path.join(os.path.join(out_dir, f"tmp{target}", smiles_file.split('/')[-1])))
            print(os.path.join(out_dir, f"tmp{target}", smiles_file.split('/')[-1]))

    aligned_dict = {'bound_pdb': [], 'smiles': []}

    for f in os.listdir(os.path.join(out_dir, f"tmp{target}")):
        if '.pdb' in f:
            aligned_dict['bound_pdb'].append(os.path.join(out_dir, f"tmp{target}", f))
            if os.path.isfile(os.path.join(out_dir, f"tmp{target}", f).replace('_bound.pdb', '_smiles.txt')):
                aligned_dict['smiles'].append(os.path.join(out_dir, f"tmp{target}", f).replace('_bound.pdb',
                                                                                               '_smiles.txt'))
            else:
                aligned_dict['smiles'].append(None)

    for aligned, smiles in list(zip(aligned_dict['bound_pdb'], aligned_dict['smiles'])):
        try:
            if smiles:
                _ = set_up(target_name=target,
                           infile=os.path.abspath(aligned),
                           out_dir=out_dir,
                           monomerize=monomerize,
                           smiles_file=os.path.abspath(smiles),
                           biomol=biomol,
                           covalent=covalent)

            else:
                _ = set_up(target_name=target,
                           infile=os.path.abspath(aligned),
                           out_dir=out_dir,
                           monomerize=monomerize,
                           biomol=biomol,
                           covalent=covalent)

        except AssertionError:
            print(aligned, "is not suitable, please consider removal or editing")
            for file in os.listdir(os.path.join(out_dir, f"tmp{target}")):
                if str(aligned) in file:
                    os.remove(os.path.join(out_dir, f"tmp{target}", str(file)))

    # Write clever method to copy in_file
    dest_dir = os.path.join(out_dir, target, 'crystallographic')
    globstr = in_file.replace('.pdb', '*')
    for file in glob.glob(globstr):
        print(file)
        shutil.copy(file, dest_dir)

    if os.path.exists(os.path.join(out_dir, f'mono{target}')):
        shutil.rmtree(os.path.join(out_dir, f'mono{target}'))

    if os.path.exists(os.path.join(out_dir, f'tmp{target}')):
        shutil.rmtree(os.path.join(out_dir, f'tmp{target}'))

    print("Files are now in a fragalysis friendly format!")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-i",
        "--in_file",
        help="Input file",
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
        "-m", "--monomerize", action="store_true", default=False, help="Monomerize input"
    )

    parser.add_argument("-t", "--target", help="Target name", required=True)

    parser.add_argument('-r', '--reference', help='Reference Structure', required=False, default=None)
    parser.add_argument('-sr', '--self_reference', action="store_true", help='Include this flag if you want each pdb file to only align to itself', required=False, default=False)
    parser.add_argument("-b", "--biomol_txt", help="Biomol Input txt file", required=False, default=None)
    parser.add_argument("-c",
                        "--covalent",
                        action="store_true",
                        help="Handle covalent bonds between ligand and target",
                        required=False,
                        default=False
                        )

    args = vars(parser.parse_args())

    in_file = args["in_file"]
    out_dir = args["out_dir"]
    monomerize = args["monomerize"]
    target = args["target"]
    biomol = args["biomol_txt"]
    covalent = args["covalent"]
    self_ref = args['self_reference']

    # Will this work?
    if self_ref:
        reference = in_file
    else if args['reference'] is None:
        reference = os.path.join(out_dir, target, 'reference.pdb')
    else:
        reference = args['reference']

    if out_dir == os.path.join("..", "..", "data", "xcimporter", "output"):
        print("Using the default input directory ", out_dir)

    if not os.path.isfile(reference):
        print(f'Cannot find file called {reference}, please make sure the path is correct (or specify another reference using -r)!')
    else:
        import_single_file(in_file=in_file,
                           out_dir=out_dir,
                           target=target,
                           monomerize=monomerize,
                           reference=reference,
                           biomol=biomol,
                           covalent=covalent,
                           self_ref=self_ref)
        print(f'File has been aligned to {reference}')


