#! /usr/env/ python

from pathlib import Path
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import pathlib


def arg_parser(args):
    parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter, 
                            description=('Takes a directory with genbank files, '
                            ' and creates a table (.tsv) with the counts of rRNAs and tRNAs.'
                            ' The 16S rRNAs must have a minimum length to be considered')
                            )
    parser.add_argument('input_dir', help = "Directory where genbank files are located")
    parser.add_argument('output_file', help="Output file, may include a path")
    parser.add_argument('--min_len', type=int, default=1400, 
                        help = "Minimum length of the 16S rRNAs to be considered")
    parser.add_argument('--glob', default='*.gbk', help = "Glob pattern to capture files")

    args = parser.parse_args()
    
    return args

def genbank_stats(gbk_file, feat_types=['CDS','rRNA', 'tRNA'], min_len_seq=1400):
    '''
        Reads a genebank file and produces a nested dictionary with stats about
        the type of features and their counts. The primary keys are the 
        features' type of interest, usually only a few.
        `min_len_seq is the min sequence length, to filter short 16S rRNAs 
        (short ones are still annotated and not marked as pseudo).
        Returns a `defaultdict` where the values of each key are class `Counter`.
    '''
    import re
    import Bio
    from Bio import SeqIO
    from Bio.Seq import Seq
    from collections import Counter, defaultdict

    parent_dict = defaultdict(list)
    counter = defaultdict()
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feat in record.features:
            for target_feat in feat_types: 
                if feat.type == target_feat:
                    for k in feat.qualifiers.keys():
                        if k == 'product' and 'pseudo' not in feat.qualifiers.keys():
                            qualifier = feat.qualifiers['product'][0]
                            if (qualifier == '16S ribosomal RNA' and 
                                len(feat.extract(record.seq)) >= min_len_seq
                               ):
                                print(record.name, qualifier,
                                      len(feat.extract(record.seq)))
                                product  = feat.qualifiers['product'][0]
                                parent_dict[feat.type].append(product)    

#                           print 16S sequences that don't pass the length filter 
                            elif qualifier == '16S ribosomal RNA':
                                print(f'{record.name}. Short 16S rRNA',
                                      qualifier, len(feat.extract(record.seq)))
                            
#                           In the case of the other annotation types, we 
#                           aren't filtering  by sequence length
                            else:
                                product  = qualifier
                                parent_dict[feat.type].append(product)
#   Create the nested dictionary with the counts of each feat type  
    for k in parent_dict.keys():
        counter[k] = Counter(parent_dict[k])
    return counter


def create_df(counter_dict, acc_name):
    '''
    Create a DataFrame from the counter dictionary coming from `genbank_stats()`.
    One DF will be created per dictionary key, where a key is the 
    feature type of interest, then DFs will be merged into one DataFrame
    '''
    import pandas as pd
    
    inner_dfs = []
    for k in counter_dict.keys():
        #Transpose it to make it "wide" 
        df = pd.DataFrame.from_dict(dict(counter_dict[k].most_common()), 
                                    orient='index',columns=[acc_name], 
                                    dtype='int').T
        inner_dfs.append(df)
        
    # concatenate along the horizontal axis. 
    df_concat = pd.concat(inner_dfs, axis=1)
    
    return df_concat

def main(args=None):

    import pandas as pd
    from pathlib import Path

    args = arg_parser(args)
    path = Path(args.input_dir)
    files = path.glob(args.glob)

    ls_dfs = []
    for gbk_file in files:
        name = gbk_file.stem
        print(name)
        counter = genbank_stats(gbk_file, feat_types=['rRNA', 'tRNA'], min_len_seq=args.min_len)
        # Some genomes may not have any tRNA or rRNA, which will produce
        # and error with create_df()
        if 'tRNA' not in counter.keys() and 'rRNA' not in counter.keys():
            print(f"{name} doesn't have any rRNA or tRNA so won't be in the final table")
            continue
        counter_df = create_df(counter, name)
        ls_dfs.append(counter_df)
        
    df_concat = pd.concat(ls_dfs, sort=True)#sort columns lexicographically
    df_concat['Total_tRNAs'] = df_concat.filter(regex='tRNA-',axis=1).sum(axis=1).astype('int')
    df_concat.rename_axis('Genome_name', inplace=True)
    df_concat.to_csv(args.output_file, sep='\t')

if __name__ == "__main__":
    main()