import unittest
from Making_Class_conversion import Ligand, set_up
import os


def test_open_file():
    """
    pdb file opens and checks that correct number of lines
    """
    new = set_up("5q1j")
    assert len(new.pdbfile) == 6779

def test_make_directory():
    """
    tests that directory to receive results has been created
    """
    new = set_up("5q1j")
    assert os.path.isdir("../results/"+str(new.pdbcode))

def test_hets_and_cons():
    """
    tests that the correct number of heteroatoms and conect lines have been found
    """
    new = set_up("5q1j")
    assert len(new.hetatms) == 326 and len(new.conects) == 24

def test_remove_nonligs():
    """
    tests that solvents and ions from crystallography have been removed
    """
    new = set_up("5q1j")
    assert len(new.final_hets) == 11

def test_wanted_ligs():
    """
    tests that the ligand identifiers at the top of the pdb file have been found
    """
    new = set_up("5q1j")
    assert len(new.wanted_ligs) != None

def test_make_mol_objs():
    """
    tests that a molecular object has been made made for each ligand
    """
    new = set_up("5qj7")
    assert len(new.mol_lst) == 4

def test_make_pdb_file():
    """
    tests that a pdb file has been made for a particular ligand with the correct number of lines (HETATM & CONECT)
    """
    set_up("5qj7")
    file = open("../results/5qj7/5qj7_JMM_A_303.pdb").readlines()
    assert len(file) == 36

def test_make_mol_file():
    """
    tests that a mol file has been made with the correct number of lines
    """
    set_up("5qj7")
    file = open("../results/5qj7/5qj7_JMM_A_303_mol.mol").readlines()
    assert len(file) == 43

def test_make_sdf_file():
    """
    test that a sdf file has been made that incorporates the different mol objects
    """
    set_up("5qj7")
    file = open("../results/5qj7/5qj7_out.sdf").readlines()
    assert len(file) == 176


if __name__ == '__main__':
    unittest.main()
