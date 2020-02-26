import urllib.request as urllib
import json
import pandas as pd

from fragalysis_api.xcglobalscripts import set_config


def xcgraphcreator(target_smiles):

    search = GraphRequest()

    # set the search smiles
    search.set_smiles_url(smiles=target_smiles)

    results = search.get_graph_json()

    results = graph_dict_to_df(results)

    return results.reset_index(drop=True)


class GraphRequest:
    def __init__(self):
        settings = set_config.ConfigSetup()

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
        with urllib.urlopen(self.smiles_url) as f:
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


def graph_dict_to_df(graph_dict):
    """
    This is the staircase to heaven
    :param graph_dict:
    :return:
    """
    a_df = pd.DataFrame()
    columns = ['type', 'insert_smiles', 'new_smiles', 'insertion']

    start = '2'

    for i1 in graph_dict[start].keys():
        for i2 in graph_dict[start][i1].keys():
            if isinstance(graph_dict[start][i1][i2], dict):
                for i3 in graph_dict[start][i1][i2].keys():
                    if isinstance(graph_dict[start][i1][i2][i3], dict):
                        for i4 in graph_dict[start][i1][i2][i3].keys():
                            if isinstance(graph_dict[start][i1][i2][i3][i4], list):
                                for i5 in graph_dict[start][i1][i2][i3][i4]:
                                    tmp_df = pd.DataFrame([[i1, i2, i5['end'], i5['change']]], columns=columns)
                                    a_df = pd.concat([a_df, tmp_df])

    return a_df
