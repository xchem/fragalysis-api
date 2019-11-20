from align import Align
from validate import Validate

if __name__ == "__main__":
    a_dir = '../anna/Hard_example/'
    out_dir = '../anna/aligned2/'

    validation = Validate(a_dir)
    validation.validate_pdbs()

    struc = Align(a_dir, pdb_ref='')
    struc.align(out_dir)