from fragalysis_api import xcgraphcreator

def xcanalyser():
    
    #def __init__(self, graphs):        
    #    self.graphs = graphs
        
    

    #base_graph = graphs
    #m = Chem.MolFromMolFile('../..data/xcimporter/output/ATAD')

    return xcgraphcreator(target_smiles='O=C(Nc1ccccc1)Nc1cccnc1')


if __name__ == '__main__':
    xcanalyser()