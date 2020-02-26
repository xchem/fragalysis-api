from fragalysis_api import xcanalyser
from rdkit import Chem


def main():

    #m = Chem.MolFromMolFile('../..data/xcimporter/output/ATAD')

    this = xcanalyser.xcgraphcreator(target_smiles='O=C(Nc1ccccc1)Nc1cccnc1')

    print(this)


if __name__ == '__main__':
    main()