from set_config import setup
import urllib
import json


class GraphRequest:
    def __init__(self):
        settings = setup()

        # get url pieces
        self.frag_url = settings.get('fragalysis', 'url')
        self.graph_url = settings.get('graph', 'search')
        self.query = settings.get('graph', 'query')

        # get full url
        self.search_url = str(self.frag_url + self.graph_url + self.query)

        # set blanks for smiles search and json to handle later
        self.smiles_url = None
        self.graph_json = None

    def set_smiles_url(self, smiles):
        # set full search url
        self.smiles_url = str(self.search_url + smiles)

    def get_graph_json(self):
        # check for a smiles url
        if not self.smiles_url:
            raise Exception('Please initiate smiles url with set_smiles_url(<smiles>)!')

        # get response from url and decode -> json
        with urllib.request.urlopen(self.smiles_url) as f:
            result = f.read().decode('utf-8')
            if not result == 'EMPTY RESULT SET':
                response = json.loads(result)

                # set json as decoded response for processing
                self.graph_json = response

        return self.graph_json

# to flatten into a list for processing
def flatten_json(y):
    out = {}

    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '_')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x
    flatten(y)

    return out


