import glob
import os
from shutil import copyfile


def to_fragalysis_dir(pdb_id, data_dir):
    """
    Creates a Fragalysis friendly directory out of a directory with the needed
    files (.sdf, .mol, _bounded.pdb and _apo.pdb). Names this new directory with the
    pdb_id name.

    :param pdb_id: ID (name) of the directory with the data
    :param data_dir: Directory which contains the data (ID directory)
    :return: Reordered directory
    """
    # Creates list of different PDBs in the directory
    pdb_aligned_files = glob.glob(os.path.join(data_dir, "*_bound.pdb"))

    # For each pdb extract the different file types of that PDB and put them in their own directory
    for pdb in [os.path.basename(pdb_file)[:-10] for pdb_file in pdb_aligned_files]:

        # Create directory if it doesn't exist
        if not os.path.exists(os.path.join(data_dir, '..', pdb_id + '_'+pdb)):
            os.makedirs(os.path.join(data_dir, '..', pdb_id + '_' + pdb))

        src = glob.glob(os.path.join(data_dir, "*" + pdb + "*"))
        new_dir = os.path.join(data_dir, '..', pdb_id + '_' + pdb)

        [copyfile(i, os.path.join(new_dir, os.path.basename(i))) for i in src]


if __name__ == '__main__':
    input_dir = '../../tests/data_for_tests/tmp/'
    to_fragalysis_dir('I_ROCK', input_dir)
