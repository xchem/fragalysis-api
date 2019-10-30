from align import Align

if __name__ == "__main__":
    dir = '/home/sabs-r3/Documents/SABS/fragalysis_api/data/ATAD2'

    struc = Align(dir)
    struc.get_files()
    struc.load()
    struc.pick_ref()
    struc.align()
    struc.split_merged()