from align import Align

if __name__ == "__main__":
    dir = '../data/ATAD2'

    struc = Align(dir, pdb_ref='')
    print(struc.get_files)
    print(struc.get_ref)
    struc.save_align()