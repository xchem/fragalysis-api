import unittest
from fragalysis_api.xcimporter import pdbquery
import os


class ValidateTest(unittest.TestCase):

    def test_query(self):
        print('testing query')
        pdb_code = '6epu'
        chain_id = 'A'
        query_obj = pdbquery.Query(pdb_code, chain_id)
        query_obj.get_matching_proteins()
        query_obj.get_ligands()
        assert(len(query_obj.match_ligs) == 41)
        query_obj.print_number_ligs()
        query_obj.view_ligands()
        query_obj.save_dictionary('anna')
        assert(os.path.exists(os.path.abspath(os.path.join('..', '..', 'data', 'xcimporter', 'other', 'anna'))))


if __name__ == '__main__':
    unittest.main()