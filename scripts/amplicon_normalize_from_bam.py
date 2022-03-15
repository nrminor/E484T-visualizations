#!/usr/bin/env python
from argparse import ArgumentParser, ArgumentTypeError


def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise ArgumentTypeError('Boolean value expected.')


def parse_process_arrays_args(parser: ArgumentParser):
    """Parses the python script arguments from bash and makes sure files/inputs are valid"""
    parser.add_argument('--bam_inpath',
                        type=str,
                        help='filepath of bamfile, should be randomized first with samtools collate or similar tool',
                        required=True)
    parser.add_argument('--bam_outpath',
                        type=str,
                        help='filepath of normalize bamfile',
                        required=True)
    parser.add_argument('--reads_per_amplicon',
                        type=int,
                        default=2000,
                        help='reads per amplicon the file should be trimmed to',
                        required=False)
    parser.add_argument('--bed_path',
                        type=str,
                        help='maximum merged read length.  Longer reads will be removed.',
                        required=True)


def get_process_arrays_args():
    """	Inputs arguments from bash
    Gets the arguments, checks requirements, returns a dictionary of arguments
    Return: args - Arguments as a dictionary
    """
    parser = ArgumentParser()
    parse_process_arrays_args(parser)
    args = parser.parse_args()
    return args


args = get_process_arrays_args()
# same arguments to a local variable by same name as the argument
bam_inpath = args.bam_inpath
bam_outpath = args.bam_outpath
reads_per_amplicon = args.reads_per_amplicon
bed_path = args.bed_path
# make sure you randomize the bam file first (i.e. samtools collate -o <outfile> <infile>)


import pandas as pd

df_bed = pd.read_csv(bed_path, sep='\t', header=None,
                     names=['SEQ_ID', 'START', 'END', 'AMPLICON_ID', 'AMPLICON_GROUP', 'SYMBOL'])
df_bed['END'] = df_bed['END']
df_bed['SIDE'] = [x.split('_')[2] for x in df_bed['AMPLICON_ID']]
df_bed['AMPLICON'] = ['_'.join(x.split('_')[:2]) for x in df_bed['AMPLICON_ID']]

group_by_list = ['AMPLICON', 'SIDE']
# because there are alternate primers, must use the max and min to encompass both sides of the primer
df_bed = df_bed.groupby(group_by_list).agg({'START': 'min',
                                            'END': 'max'}).reset_index()
# pivot so the start and ends of the forward and reverse primers are in the same row,
df_bed_pivot = df_bed.pivot_table(index='AMPLICON', columns='SIDE', values=['START', 'END']).reset_index()
# rename columns ass needed
df_bed_pivot.columns = df_bed_pivot.columns.to_series().str.join('_')
df_bed_pivot.rename(columns={'AMPLICON_': 'AMPLICON'}, inplace=True)
df_bed_pivot.sort_values(by=['START_LEFT'], inplace=True)
# df_bed_pivot.reset_index(inplace=True)
# df_bed_pivot.drop(columns=['index'], inplace = True)
prior_start = list(df_bed_pivot['START_LEFT'])
prior_start.insert(0, -100000)
prior_start.pop()
next_start = list(df_bed_pivot['START_LEFT'])
next_start.pop(0)
max_next = max(next_start) * 1000
next_start.append(max_next)
df_bed_pivot['PRIOR_START'] = prior_start
df_bed_pivot['NEXT_START'] = next_start
df_bed_pivot = df_bed_pivot[['AMPLICON', 'START_LEFT', 'PRIOR_START', 'NEXT_START']]
bed_dict = {}
amplicon_count = {}
for index, row in df_bed_pivot.iterrows():
    bed_dict[row['AMPLICON']] = [row['PRIOR_START'], row['NEXT_START']]
    amplicon_count[row['AMPLICON']] = 0

import pysam

bf_in = pysam.AlignmentFile(bam_inpath, "rb")
bf_out = pysam.AlignmentFile(bam_outpath, "wb", template=bf_in)
i = 0
for r in bf_in.fetch(until_eof=True):
    i = i + 1
    if i % 20000 == 0:
        print('PROCESSED READS = {0}, Start Position = {1}'.format(i, start))
    mapped_read = str(r).split('\t')
    start = int(mapped_read[3])
    # end = mapped_read[8]
    for key, value in bed_dict.items():
        if value[0] < start < value[1]:
            amplicon_count[key] = amplicon_count[key] + 1
            if amplicon_count[key] > reads_per_amplicon:
                break
            bf_out.write(r)
            break
bf_out.close()
bf_in.close()
