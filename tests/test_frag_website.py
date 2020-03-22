import unittest
from urllib.request import urlopen
from urllib.error import URLError

from fragalysis_api.xcglobalscripts import set_config
from fragalysis_api.xcextracter import frag_web_live


class WebFragalysis(unittest.TestCase):

    def test_fragalysis_connection(self):

        def frag_connection():
            settings = set_config.ConfigSetup()
            url = settings.get('fragalysis', 'url')

            try:
                urlopen(url, timeout=1)
                status = True
            except URLError as err:
                status = False

            return status

        if frag_web_live.can_connect():
            print('Can connect to Fragalysis')
        else:
            print('Cannot connect to Fragalysis')
        
        self.assertTrue(frag_connection())

    def test_can_connect_function(self):
            
        self.assertIsNotNone(frag_web_live.can_connect())


if __name__ == '__main__':
    unittest.main()
