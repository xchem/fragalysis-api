from align import Align
from validate import Validate

if __name__ == "__main__":
    a_dir = '../data/ATAD2/'

    #a_dir = '../data/Hard_example/'


    validation = Validate(a_dir)
    validation.validate_pdbs

    struc = Align(a_dir, pdb_ref='')
    struc.align()