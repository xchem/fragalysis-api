import unittest
from urllib.request import urlopen
from urllib.error import URLError
from fragalysis_api.xcglobalscripts import set_config


class test_fragalysis(unittest.TestCase):

    def test_fragalysis_connection(self):

        def can_connect():
            settings = set_config.ConfigSetup()
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
        
        self.assertTrue(can_connect())


if __name__ == '__main__':
    unittest.main()
