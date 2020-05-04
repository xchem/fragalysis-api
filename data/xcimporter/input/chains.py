import os
import json

non_ligs = json.load(open('non_ligs.json', 'r'))

def get_lig_chains(file):
    pdb = open(os.path.join('Mpro_PDB', file)).readlines()
    ligs = {}
    for line in pdb:
        if line.startswith('HET '):
            lig, chain =  line[7:10], line[12:13]
            if lig not in non_ligs:
                if lig not in ligs:
                    ligs[lig] = [chain]
                else:
                    ligs[lig].append(chain)
    return ligs


def ask_which_lig(ligs):
    all_ligs = [i for i in ligs]
    print('Many ligs found in', file)
    print(set(all_ligs))
    lig = input('Which ligand would you like to view in fragalysis?')
    chains_to_keep = [ligs[lig]]
    if input('Would you like to view another ligand in fragalysis?') == 'y':
        lig2 = input('Which other ligand would you like to view in fragalysis?')
        chains_to_keep.append(ligs[lig2])


for file in os.listdir('Mpro_PDB'):
    ligs = get_lig_chains(file)
    print(ligs)
    # try:
    #     if len(ligs) > 1:
    #         chains = ask_which_lig(ligs)
    #     else:
    #         chains = [ligs[i] for i in ligs]
    # except IndexError:
    #     print(file)