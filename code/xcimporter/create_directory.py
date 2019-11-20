import glob
import os
from shutil import copyfile


def create_directory(PDBcode, data_dir):
    pdb_aligned_files = glob.glob(os.path.join(data_dir, "*_bound.pdb"))

    for pdb in [os.path.basename(pdb_file)[:-10] for pdb_file in pdb_aligned_files]:
        if not os.path.exists(os.path.join(data_dir, '..', PDBcode+'_'+pdb)):
            os.makedirs(os.path.join(data_dir, '..', PDBcode+'_'+pdb))

        src = glob.glob(os.path.join(data_dir, "*"+pdb+"*"))

        [copyfile(i, os.path.join(data_dir, '..', PDBcode+'_'+pdb, os.path.splitext(os.path.basename(i))[0])) for i in src]


if __name__ == '__main__':
    input_dir = '../../data/xcimporter/output/ATAD2/tmp/'
    create_directory('I_ROCK', input_dir)
