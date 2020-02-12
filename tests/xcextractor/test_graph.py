import unittest
from graph import GraphRequest

class GraphRequestTest(unittest.TestCase):

    def test_GraphRequest_exists(self):
        g = GraphRequest()
        self.assertIsNotNone(g)

    def test_set_smiles_url(self):
        g = GraphRequest()
        g.set_smiles_url(smiles = 'O=C(Nc1ccccc1)Nc1cccnc1')
        self.assertIsNotNone(g.smiles_url)

    def test_get_graph_json(self):
        g = GraphRequest()
        g.set_smiles_url(smiles = 'O=C(Nc1ccccc1)Nc1cccnc1')
        result = g.get_graph_json()
        self.assertIsNotNone(result)

if __name__ == '__main__':
    unittest.main()