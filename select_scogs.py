#! /usr/bin/env python

from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import shutil

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import textwrap

def arg_parser(args):
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, 
                            description='Selects single-copy orthogroups present in all the'
                                         ' genomes of a given node of the'
                                         ' species tree created by Orthofinder')

    parser.add_argument('-hogs_dir', default="Phylogenetic_Hierarchical_Orthogroups",
                        help = 'Directory containing all the Nx.tsv files.')
    parser.add_argument('-node_to_use', 
                        help = 'Nx.tsv file to be used', default="N0.tsv")
    parser.add_argument('-ogs_path', default="Orthogroup_Sequences",
                        help = 'Orthogroups sequences path')
    parser.add_argument('-out_dir', default='Single_copy_OGs',
                        help = 'Directory where the scogs will be placed')
    args = parser.parse_args()
    
    return args

def select_scogs(hogs_dir, node_to_use):
    '''
    Takes a specific 'Nx.tsv' file and selects the single copy orthogroups present
    in all species (core SCOGs).
    '''

    hogs_dir = Path(hogs_dir)
    n_df = pd.read_csv(hogs_dir.joinpath(node_to_use), sep='\t', na_filter=False)

    og_dict = defaultdict()
    genome_dict = defaultdict(list)
    for i, row in n_df.iterrows():
        genome_dict = defaultdict(list)
        for asm, gene in row[3:].items():
            if gene != '':
                genome_dict[asm].append(gene.split(','))
        og_dict[row['OG']] = genome_dict
    # select OGs present in all the genomes
    all_asms = np.array(n_df.columns[3:].to_list())
    og_present_inall = []
    i = 0
    for k,v in og_dict.items():
        i += 1
        if np.isin(all_asms, np.array(list(v.keys()))).all():
            og_present_inall.append(k)
    print("Number of OGs present in all species: ", len(og_present_inall))
    
    # Define core single copy OGs
    single_copy_core_ogs = []
    for og in og_present_inall:
        swt = False
        for k,v in og_dict[og].items():
    #         unnest the nested list
            genes = [gene for ls in v for gene in ls]
            if len(genes) == 1:
                swt = True
            if len(genes) > 1:
                swt = False
                break
        if swt:
            single_copy_core_ogs.append(og)
            
    print("Number of core single-copy OGs: ",len(single_copy_core_ogs))
    return(single_copy_core_ogs)


def copy_scogs(hogs_dir, node_to_use, ogs_path,
               out_dir):
    '''
    Takes a list with the single-copy OGs and copies them to 'outdir'
    '''
    
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    
    scog_list = select_scogs(hogs_dir, node_to_use)
    scogs_fas = (pd.Series(scog_list) + ".fa").to_list()
    ogs_path = Path(ogs_path)
    i = 0
    for fa in ogs_path.glob("*.fa"):
        if fa.name in scogs_fas:
            i += 1
            shutil.copy(fa, out_dir)

    print(f'{i} fa files were copied to {out_dir}')
    
    assert i == len(scog_list)
    
    return 0

def main(args=None):
    args = arg_parser(args)
#     print(args)
    copy_scogs(args.hogs_dir, args.node_to_use, args.ogs_path, args.out_dir)

if __name__ == "__main__":
    main()