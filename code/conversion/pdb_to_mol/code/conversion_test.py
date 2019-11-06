import unittest
from Making_Class_conversion import Ligand, set_up
import os


def test_open_file():
    new = set_up("5q1j")
    assert len(new.pdbfile) == 6779

def test_make_directory():
    new = set_up("5q1j")
    assert os.path.isdir("../results/"+str(new.pdbcode))

def test_hets_and_cons():
    new = set_up("5q1j")
    assert len(new.hetatms) == 326 and len(new.conects) == 24

def test_remove_nonligs():
    new = set_up("5q1j")
    assert len(new.final_hets) == 11

def test_wanted_ligs():
    new = set_up("5q1j")
    assert len(new.wanted_ligs) != None

def test_make_mol_objs():
    new = set_up("5qj7")
    assert len(new.mol_lst) == 4

def test_make_pdb_file():
    set_up("5qj7")
    file = open("../results/5qj7/5qj7_JMM_A_303.pdb").readlines()
    assert len(file) == 36

def test_make_mol_file():
    set_up("5qj7")
    file = open("../results/5qj7/5qj7_JMM_A_303_mol.mol").readlines()
    assert len(file) == 43

def test_make_sdf_file():
    set_up("5qj7")
    file = open("../results/5qj7/5qj7_out.sdf").readlines()
    assert len(file) == 176


if __name__ == '__main__':
    unittest.main()
