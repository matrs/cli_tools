#! /usr/bin/env python


from pathlib import Path
import subprocess as sp
from collections import defaultdict
import re, os
from ete3 import Tree
import pickle
import random
import string

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def arg_parser(args):
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, 
                            description='Takes a gene tree and the species tree from Orthofinder'
                                         ' and runs Ranger-DTL and Aggregate Ranger.  All the leaves '
                                         ' in the species and gene trees get transformed to random strings'
                                         ' because Ranger-DTL issues. The aggregated files get converted' 
                                         ' back and a summary table gets created. A species Tree and a '
                                         ' pickled dictionary are hardcoded')

    parser.add_argument('species_tree', help = 'Previously edited species tree from Orthofinder')
    parser.add_argument('gene_tree', help = 'Gene tree from Orthofinder.')
    parser.add_argument('snames_dict', 
                                    help='Species names to random string dictionary (pickle).'
                                    ' ')
    parser.add_argument('out_dir', help = 'Output directory where the tree and dictionary will be placed.')
    args = parser.parse_args()
    
    return args


def create_gene_dict(tree_file, names_dict, out_dir):
    '''
    Takes as input a gene tree, not a species tree
    '''
    ete_tree = Tree(tree_file, format= 1)
    sname_gene_to_rand = defaultdict()
    rand_to_sname_gene = defaultdict()
    for node in ete_tree.traverse():
        if node.is_leaf():
            for k,v in names_dict.items():
                if re.match(k, node.name):
                    gene_name = '_'.join(node.name.split('_')[-2:])
                    sname_gene_to_rand[f'{k}_{gene_name}'] = f'{v}_{gene_name}' 
                    rand_to_sname_gene[f'{v}_{gene_name}' ] = f'{k}_{gene_name}'
    # Create a mapping for the gene tree: random ==> spname_gene_name
    pickle_out = out_dir.joinpath(f'new_{Path(tree_file).stem}.pkl')
    with open(pickle_out, 'wb') as fh:
        pickle.dump(rand_to_sname_gene, fh)

    return(sname_gene_to_rand)


def convert_names(tree_file, names_dict, out_dir):
    '''Takes as input a gene tree, not a species tree'''
    ete_tree = Tree(tree_file, format= 1)
    with open(names_dict, 'rb') as handle:
        names_dict = pickle.load(handle)
    sname_gene = create_gene_dict(tree_file, names_dict, out_dir)
#     print(sname_gene)
    for node in ete_tree.traverse():
        if node.is_leaf():
            node.name = sname_gene[node.name]
    return(ete_tree)


ranger = "/home/jmaturana/Software/RangerDTL_Linux/CorePrograms/Ranger-DTL.linux"
aggregate = "/home/jmaturana/Software/RangerDTL_Linux/CorePrograms/AggregateRanger.linux"
D,T,L = [2,3,2],[3,3,4],[1,1,1]
n_seeds = 10
print(f'Number of different seeds {n_seeds}')


def main(args=None):
    args = arg_parser(args)
    stree = args.species_tree
    og_tree = args.gene_tree
    out_dir = Path(args.out_dir)
    agg_dir = out_dir.joinpath('AggregateRanger')
    agg_dir.mkdir(parents=True, exist_ok=True)
    # Create the output dir for each OG
    stem_og = Path(og_tree).stem
    tree_out = out_dir.joinpath(stem_og) # This needs to be Path to use .mkdir
    # new gene_tree to concatenate it with the species tree
    gtree = convert_names(og_tree, args.snames_dict, out_dir)
    gtree_outf = out_dir.joinpath(f'new_{stem_og}.txt')
    gtree.write(format=1, outfile=gtree_outf.as_posix(), dist_formatter='%f')
    input_tree = out_dir.joinpath(f"input_{stem_og}.nwk")
    output_dir = out_dir.joinpath(stem_og)
    with open(input_tree, 'w') as fh:
    #Add the "root" node, n0, given that ete takes it off regardless the format
    # N0 was added the same way to the species tree
        print(Path(stree).read_text(), Path(gtree_outf).read_text().replace(';','n0;',1), 
                file=fh, sep='')
    #Remove the edited gene tree after concatenation
    os.remove(gtree_outf)
    print(stem_og)
    i = 0
    for seed in range(1, n_seeds+1):
        # print('Looping, seed = ', seed)
        for d,t,l in zip(D,T,L):
            tree_out.mkdir(exist_ok=True)
            # out_dir.joinpath('log').mkdir(exist_ok=True)
            # log_out = out_dir.joinpath('log').as_posix()
            i += 1
            ranger_cmd = (f"{ranger} --seed {seed} -i {input_tree} -D {d} -T {t} -L {l} "
                            f"-o {output_dir.as_posix()}/out{i}")
            sp.check_output(ranger_cmd, shell=True, universal_newlines=True)
        ## aggregate over all the weights and seeds ##
    agg_outf = agg_dir.joinpath(f"agg_{stem_og}.txt")
    agg_cmd = f"{aggregate} {output_dir.as_posix()}/out > {agg_outf.as_posix()}"
    sp.check_output(agg_cmd, shell=True, universal_newlines=True)
    # remove the concatenated tree
    os.remove(input_tree)


if __name__ == "__main__":
    main()