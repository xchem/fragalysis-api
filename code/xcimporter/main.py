from align import Align
from validate import Validate
import os

if __name__ == "__main__":
    a_dir = '../../data/anna/input/ATAD2/'
    out_dir = '../../data/output/ATAD2/'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    validation = Validate(a_dir)
    validation.validate_pdbs()

    struc = Align(a_dir, pdb_ref='')
    struc.align(out_dir)