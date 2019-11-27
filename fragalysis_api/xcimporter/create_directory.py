import glob
import os
from shutil import copyfile


def create_directory(PDBcode, data_dir):
    """
    Creates a fragalysis friendly directory out of a directory with the needed
    files (.sdf, .mol, _bounded.pdb and _apo.pdb). Names this new directory with the
    PDBcode name.

    :param PDBcode: ID (name) of the directory with the data
    :param data_dir: Directory which contains the data (ID directory)
    :return:
    """
    # Creates list of different pdbs in the directory
    pdb_aligned_files = glob.glob(os.path.join(data_dir, "*_bound.pdb"))

    # For each pdb extract the different file types of that pdb and put them in their own directory
    for pdb in [os.path.basename(pdb_file)[:-10] for pdb_file in pdb_aligned_files]:
        # Create directory if it doesn't exist
        if not os.path.exists(os.path.join(data_dir, '..', PDBcode+'_'+pdb)):
            os.makedirs(os.path.join(data_dir, '..', PDBcode+'_'+pdb))


        src = glob.glob(os.path.join(data_dir, "*"+pdb+"*"))
        new_dir = os.path.join(data_dir, '..', PDBcode+'_'+pdb)

        [copyfile(i, os.path.join(new_dir, os.path.basename(i))) for i in src]


if __name__ == '__main__':
    input_dir = '../../data/xcimporter/output/ATAD2/tmp/'
    create_directory('I_ROCK', input_dir)
