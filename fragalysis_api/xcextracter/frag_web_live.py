from urllib.request import urlopen
from urllib.error import URLError
from fragalysis_api import ConfigSetup


def can_connect():
    '''
    Function to check if the fragalysis website specified in the config.ini is live
    return status: boolean specifying if fragalysis is live
    '''

    settings = ConfigSetup()
    #defining the fragalysis url
    url = settings.get('fragalysis', 'url')

    status = ''
    try:
        #checking it responses in a certain time
        urlopen(url, timeout=1)
        status = True
    except URLError as err:
        status = False

    return status