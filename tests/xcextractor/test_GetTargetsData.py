import unittest
from fragalysis_api.xcextracter import getdata


class GetTargetsDataTest(unittest.TestCase):

    def test_GetTargetsData_exists(self):
        search = getdata.GetTargetsData()
        self.assertIsNotNone(search)
    
    def test_set_target_name_url(self):
        search = getdata.GetTargetsData()
        search.set_target_name_url('ATAD')
        self.assertIsNotNone(search.target_name_url)
    
    def test_get_target_json(self):
        search = getdata.GetTargetsData()
        search.set_target_name_url('ATAD')
        self.assertIsNotNone(search.get_target_json)

    def test_get_target_id_list(self):
        search = getdata.GetTargetsData()
        search.set_target_name_url('ATAD')
        self.assertIsNotNone(search.get_target_id_list)


if __name__ == '__main__':
    unittest.main()
