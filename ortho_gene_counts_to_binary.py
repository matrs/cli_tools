#! /usr/bin/env python

#python ortho_gene_counts_to_binary -h

from pathlib import Path
import pandas as pd
import numpy as np
from argparse import ArgumentParser

def genecounts_to_binary(gene_count_tsv, out_file):
    '''
    Takes in the `Orthogroups.GeneCount.tsv` and converts it into
    a binary table where 1s indicate presence of an orthogroup.
    Also, takes out the column 'Total' because isn't needed.
    Output: out_file, name of the file to be written. It May include the path
    '''
    
    # lambda col: col != 'Total' == > function will be called against the columns' names. 
    # This excludes the column 'Total' because makes `count` throw an error and it's not 
    # needed anyways.
    input = Path(gene_count_tsv)
    output = Path(out_file)
    orthogroups_df = pd.read_csv(input, index_col= 0, sep='\t', 
                                 usecols = lambda col: col != 'Total')
    orth_df_bin = (orthogroups_df >= 1).astype(np.int8)
    orth_df_bin.to_csv(output, sep='\t', index_label='Orthogroup')
    
    return(print(f'Binary table written into {out_file}'))

def arg_parser(args):
    parser = ArgumentParser(prog='Orthofinder Gene counts to binary', 
                            description=('Creates a binary table where 1s indicate '
                                         'presence of an specific orthogroup'))
    parser.add_argument('gene_counts',
                        help = ('Orthogroups.GeneCount.tsv generated by Orthofinder, '
                        'Gene count table for each orthogroup coming from Orthofinder')
                        )
    parser.add_argument('out_file', 
                        help = 'Name of the file to be written, may include the path')
    args = parser.parse_args()
    
    return args

def main(args=None):
    args = arg_parser(args)
#     print(args)
    genecounts_to_binary(args.gene_counts, args.out_file)

if __name__ == "__main__":
    main()
