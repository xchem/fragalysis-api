import unittest
from fragalysis_api import GetMoleculesData

class GetMoleculesDataTest(unittest.TestCase):

    def test_GetMoleculesData_exists(self):
        search = GetMoleculesData()
        self.assertIsNotNone(search)
    
    def test_get_target_id(self):
        search = GetMoleculesData()
        search.set_target_id('ATAD')
        self.assertIsNotNone(search.get_target_id)
    
    def test_set_molecule_url(self):
        search = GetMoleculesData()
        search.set_target_id('ATAD')
        search.set_molecule_url()
        self.assertIsNotNone(search.get_molecule_url)
    
    def test_set_mol_data(self):
        search = GetMoleculesData()
        search.set_target_id('ATAD')
        search.set_molecule_url()
        search.set_mol_data()
        self.assertIsNotNone(search.get_mol_data)

    def test_set_complete_mol_data(self):
        search = GetMoleculesData()
        search.set_target_id('ATAD')
        search.set_molecule_url()
        search.set_mol_data()
        search.set_complete_mol_data()
        self.assertIsNotNone(search.get_complete_mol_data)


if __name__ == '__main__':
    unittest.main()