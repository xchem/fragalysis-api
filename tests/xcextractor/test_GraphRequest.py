import unittest
from fragalysis_api import GraphRequest


class test_graphRequest(unittest.TestCase):

    def test_graphRequest(self):
        search = GraphRequest()
        self.assertIsNotNone(search)

    def test_set_smiles_url(self):
        target_smiley='O=C(Nc1ccccc1)Nc1cccnc1'
        search = GraphRequest()
        search.set_smiles_url(smiles=target_smiley)
        self.assertIsNotNone(search.smiles_url)

    def test_get_graph_json(self):
        target_smiley='O=C(Nc1ccccc1)Nc1cccnc1'
        search = GraphRequest()
        search.set_smiles_url(smiles=target_smiley)
        self.assertIsNotNone(search.get_graph_json)

if __name__ == '__main__':
    unittest.main()