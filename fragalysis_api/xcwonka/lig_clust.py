import os
import argparse

from fragalysis_api.xcwonka import funcs_for_lig_clust

def lig_cluster(DATA_DIRECTORY: str = os.path.abspath('data')):
    '''
    clusters by ligands
    :param DATA_DIRECTORY: the directory in which the data director is in
    :return output: clusters of ligands 
    '''
   
    print('directory: ', DATA_DIRECTORY)

    for dir in os.listdir(DATA_DIRECTORY):

        print('directory: ', dir)
        mols = []
        identifiers = []

        for file in os.listdir(os.path.join(DATA_DIRECTORY, dir)):

            if 'bound' in file:

                print('bound_file: ', file)
                mols.append(funcs_for_lig_clust._parse_pdb(os.path.join(DATA_DIRECTORY, dir, file)))
                identifiers.append(str(file))

    output = funcs_for_lig_clust.run_lig_cluster(mols, identifiers)

    return output

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-id', '--user_id', required=True,
                        help='Description for foo argument')
    args = vars(parser.parse_args())

    user_id = args['user_id']

    in_dir =  os.path.join('..', '..', 'data', 'xcimporter', 'output', user_id)

    trial = lig_cluster(in_dir)

    print('final output: ', trial)
    print('The ligands have been clustered!!!!')