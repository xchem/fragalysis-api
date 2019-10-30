import glob
import Bio.PDB
from pymol import cmd

lst = glob.glob("*.pdb")
lst.sort()

parser = Bio.PDB.PDBParser()


class PDBFile:
    def __init__(self, name, length, resolution):
        self.name = name
        self.length = length
        self.resolution = resolution


for i in glob.glob(".pdb"):
    cmd.load(i)
    cmd.align(lst[i],lst[0])

cmd.multisave("5qj_all_aligned.pdb")