import argparse
import glob
import os
import shutil

from fragalysis_api import Align, set_up, convert_small_AA_chains, copy_extra_files, Sites, contextualize_crystal_ligands


def import_single_file(in_file, out_dir, target, reduce_reference_frame, reference_pdb, biomol=None, covalent=False, self_ref=False, max_lig_len=0):
    '''Formats a PDB file into fragalysis friendly format.
    1. Validates the naming of the pdbs.
    2. It aligns the pdbs (_bound.pdb file).
    3. It cleans the pdbs (removes solvents and ions).
    4. Splits the pdbs into ligands (.mol, .sdf and .pdb formats) and protein files (_apo.pdb file).
    5. Orders all the files into the correct directory structure required for fragalysis.
    :param in_file: Filepath of pdb you with to align (where additional files are stored in same directory)
    :param out_dir: Directory containing processed pdbs (will be created if it doesn't exists)
    :param target: Name of the folder to be created inside out_dir
    :param reduce_reference_frame: Bool, if True, will attempt to split pdb files into seperate chains
    :param reference_pdb: Name of the Reference pdb to align in_file to. If called from command line this will default to out_dir/target/reference.pdb
    :param biomol: plain-text file containing header information about the bio-molecular
        context of the pdb structures. If provided the contents will be appended to the top of the _apo.pdb files
    :param covalent: Bool, if True, will attempt to convert output .mol files to account for potential covalent attachments
    :param self_ref: Bool, if True, the import single file will align to itself.
    :max_lig_len: Integer, If >0 will convert all chains with fewer than max_lig_len residues to HETATM with the name LIG. [Currently broken, yikes]
    :return: Hopefully, beautifully aligned files that be used with the fragalysis loader :)
    '''

    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    if not os.path.isdir(os.path.join(out_dir, f"tmp{target}/")):
        os.makedirs(os.path.join(out_dir, f"tmp{target}/"))

    in_dir = os.path.dirname(in_file)

    pdb_smiles_dict = {'pdb': [], 'smiles': []}

    # Experimental - if option is used then chains are converted.
    if int(max_lig_len) > int(0):
        print(
            f'EXPERIMENTAL: Converting all chains with less than {max_lig_len} residues to HETATM LIG')
        out = os.path.join(out_dir, f'maxliglen{target}/')
        if not os.path.isdir(out):
            os.makedirs(out)
        convert_small_AA_chains(in_file=in_file, out_file=os.path.join(
            out, in_file), max_len=max_lig_len)
        copy_extra_files(in_file=in_file, out_dir=out)
        in_dir = out

    pdb_smiles_dict['pdb'].append(in_file)
    if os.path.isfile(in_file.replace('.pdb', '_smiles.txt')):
        pdb_smiles_dict['smiles'].append(
            in_file.replace('.pdb', '_smiles.txt'))
    else:
        pdb_smiles_dict['smiles'].append(None)

    print(f"Aligning to Reference: {reference_pdb}")
    structure = Align(in_dir, "", rrf=reduce_reference_frame,
                      refset=False)

    for i in pdb_smiles_dict['pdb']:
        if self_ref:
            structure.align_to_reference(
                i, i, out_dir=os.path.join(out_dir, f"tmp{target}"))
        else:
            structure.align_to_reference(
                i, reference_pdb=reference_pdb, out_dir=os.path.join(out_dir, f"tmp{target}"))

    for smiles_file in pdb_smiles_dict['smiles']:
        if smiles_file:
            shutil.copyfile(smiles_file, os.path.join(os.path.join(
                out_dir, f"tmp{target}", smiles_file.split('/')[-1])))
            print(os.path.join(
                out_dir, f"tmp{target}", smiles_file.split('/')[-1]))

    aligned_dict = {'bound_pdb': [], 'smiles': []}
    for f in os.listdir(os.path.join(out_dir, f"tmp{target}")):
        if '.pdb' in f:
            aligned_dict['bound_pdb'].append(
                os.path.join(out_dir, f"tmp{target}", f))
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
                           rrf=reduce_reference_frame,
                           smiles_file=os.path.abspath(smiles),
                           biomol=biomol,
                           covalent=covalent,
                           keep_headers=True)

            else:
                _ = set_up(target_name=target,
                           infile=os.path.abspath(aligned),
                           out_dir=out_dir,
                           rrf=reduce_reference_frame,
                           biomol=biomol,
                           covalent=covalent,
                           keep_headers=True)

        except AssertionError:
            print(aligned, "is not suitable, please consider removal or editing")
            for file in os.listdir(os.path.join(out_dir, f"tmp{target}")):
                if str(aligned) in file:
                    os.remove(os.path.join(out_dir, f"tmp{target}", str(file)))

    # Write clever method to copy in_file
    dest_dir = os.path.join(out_dir, target, 'crystallographic')
    globstr = in_file.replace('.pdb', '*')
    for file in glob.glob(globstr):
        shutil.copy(file, dest_dir)

    # Time to use a for loop?
    clean_up = [os.path.join(out_dir, f'maxliglen{target}'), os.path.join(
        out_dir, f'mono{target}'), os.path.join(out_dir, f'tmp{target}')]
    [shutil.rmtree(x) for x in clean_up if os.path.exists(x)]

    # Cut Maps
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
        "-rrf", "--reduce_reference_frame", action="store_true", default=False, help="Reduce Reference Frame"
    )

    parser.add_argument("-t", "--target", help="Target name", required=True)

    parser.add_argument('-r', '--reference_pdb',
                        help='Reference PDB Structure', required=False, default=None)
    parser.add_argument('-sr', '--self_reference', action="store_true",
                        help='Include this flag if you want each pdb file to only align to itself', required=False, default=False)
    parser.add_argument(
        "-b", "--biomol_txt", help="Biomol Input txt file", required=False, default=None)
    parser.add_argument("-c",
                        "--covalent",
                        action="store_true",
                        help="Handle covalent bonds between ligand and target",
                        required=False,
                        default=False
                        )
    parser.add_argument("-mll",
                        "--max_lig_len",
                        help="Int, Convert all chains shorter than max_lig_len to HETATM LIG",
                        required=False,
                        default=0)

    parser.add_argument(
        "-cs",
        "--cluster_sites",
        action="store_true",
        help='Include this flag if you would like to automatically assign center of mass sites as site labels',
        required=False,
        default=False
    )

    parser.add_argument(
        "-cs_com",
        "--cluster_sites_com",
        help="Tolerance value for creating new clusters for centre of mass sites",
        type=float,
        default=5.00
    )

    parser.add_argument(
        "-cs_other",
        "--cluster_sites_other",
        help="Tolerance value for creating new clusters for non centre of mass sites",
        type=float,
        default=1.00
    )

    args = vars(parser.parse_args())

    in_file = args["in_file"]
    out_dir = args["out_dir"]
    reduce_reference_frame = args["reduce_reference_frame"]
    target = args["target"]
    biomol = args["biomol_txt"]
    covalent = args["covalent"]
    self_ref = args['self_reference']
    mll = args['max_lig_len']
    cs = args['cluster_sites']
    cs_com = args['cluster_sites_com']
    cs_other = args['cluster_sites_other']

    # Will this work?
    if self_ref:
        reference_pdb = in_file
    elif args['reference_pdb'] is None:
        reference_pdb = os.path.join(out_dir, target, 'reference.pdb')
    else:
        reference_pdb = args['reference_pdb']

    if out_dir == os.path.join("..", "..", "data", "xcimporter", "output"):
        print("Using the default input directory ", out_dir)

    if not os.path.isfile(reference_pdb):
        print(
            f'Cannot find file called {reference_pdb}, please make sure the path is correct (or specify another reference using -r)!')
    else:
        import_single_file(in_file=in_file,
                           out_dir=out_dir,
                           target=target,
                           reduce_reference_frame=reduce_reference_frame,
                           reference_pdb=reference_pdb,
                           biomol=biomol,
                           covalent=covalent,
                           self_ref=self_ref,
                           max_lig_len=mll)
        print(f'File has been aligned to {reference_pdb}')
        if cs:
            folder = os.path.join(out_dir, target)
            site_obj = Sites.from_folder(folder, recalculate=False)
            site_obj.cluster_missing_mols(
                com_tolerance=cs_com, other_tolerance=cs_other)
            site_obj.to_json()
            contextualize_crystal_ligands(folder=folder)
            site_obj.apply_to_metadata()
