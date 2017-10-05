#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 16:20 07/08/2017 2017 

"""
import sys
import os
from subprocess import check_call

import pandas as pd

__author__ = 'tbrittoborges'
path = os.path


def read_gff(filename):
    names = 'seqname source feature start end score strand frame attribute'.split()
    dtypes = {k: v for k, v in zip(names,
                                   'str str str int int float str str str'.split())}

    return pd.read_table(
        filename,
        comment='#',
        dtype=dtypes,
        names=names,
        na_values='.'
    )


def create_file(filename, command):
    command = command.split()
    with open(filename, 'w') as open_f:
        check_call(
            command,
            stdout=open_f
        )


def split_features(features=None):
    if features is None:
        features = ['gene', 'CDS', 'three_prime', 'five_prime', 'exon']

    for feature in features:
        temp = df.query('features == "{}"'.format(feature))
        temp.to_table(
            "{}/{}.gff".format(assembly, feature),
            sep='\t',
            na_values='.',
            index=False,
            header=False
        )


assembly, extension = path.splitext(path.basename(sys.argv[1]))
support_ext = ['.gff', '.gff3']
if extension not in support_ext:
    # no file validation here, only naively checks the file extension
    error = 'Unsupported file extension: {} Supported extensions are: {}'
    raise RuntimeError(error.format(extension, ", ".join(support_ext)))

try:
    os.makedirs(assembly)
except OSError as e:
    pass  # dir exists

df = read_gff(sys.argv[1])
split_features()

# filenames = ["{}/{}.gff".format(assembly, 'intron'),
#              "{0}/{0}.'genome'".format(assembly),
#              "{}/intergenic.bed".format(assembly)
#              ]
#
# cmds = ['bedtools subtract -s -a genes.gff -b exon.gff |sed -e "s/\tgene\t/\tintron\t/"',
#         'mysql - -user = genome - -host = genome - mysql.cse.ucsc.edu - A - e "select chrom, '
#         'size from ${ASSEMBLY}.chromInfo',
#         'bedtools complement -i genes.gff -g {}.genome'
#         ]
#
# for command, filename in zip(cmds, filenames):
#     create_file(filename, command)

with open("{}/{}.gff".format(assembly, 'intron'), 'w') as open_f:
    check_call(
        'bedtools subtract -s -a genes.gff -b exon.gff |sed -e "s/\tgene\t/\tintron\t/"'.split(),
        stdout=open_f
    )

with open("{0}/{0}.'genome'".format(assembly), 'w') as open_f:
    check_call(
        'mysql - -user = genome - -host = genome - mysql.cse.ucsc.edu - A - e "select chrom, '
        'size from ${ASSEMBLY}.chromInfo"'.split(),
        stdout=open_f
    )

with open("{}/intergenic.bed".format(assembly), 'w') as open_f:
    check_call(
        'bedtools complement -i genes.gff -g {}.genome'.format(assembly).split(),
        stdout=open_f
    )
