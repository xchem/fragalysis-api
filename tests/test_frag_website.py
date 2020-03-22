import unittest
from urllib.request import urlopen
from urllib.error import URLError
from fragalysis_api import ConfigSetup
from fragalysis_api import can_connect


class test_fragalysis(unittest.TestCase):

    def test_fragalysis_connection(self):

        def frag_connection():
            settings = ConfigSetup()
            url = settings.get('fragalysis', 'url')
            status = ''
            try:
                urlopen(url, timeout=1)
                status = True
            except URLError as err:
                status = False

            return status
        if can_connect():
            print('Can connect to Fragalysis')
        else:
            print('Cannot connect to Fragalysis')
        
        self.assertTrue(frag_connection())

    def test_can_connect_function(self):
            
        self.assertIsNotNone(can_connect())


if __name__ == '__main__':
    unittest.main()