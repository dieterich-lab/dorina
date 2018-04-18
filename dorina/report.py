#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 13:48 16/04/2018 2018 

"""
import multiprocessing

import numpy as np
import pandas as pd
import pybedtools
from bokeh.models import ColumnDataSource, LabelSet
from bokeh.plotting import figure, show


def load_gff_to_df(filename, filetype='bigNarrowPeak'):
    """
    Loads a clip data like file to pandas.Dataframe


    :param filename:
    :param filetype:
    :return:

    ..note::
        see https://genome.ucsc.edu/goldenPath/help/bigNarrowPeak.html
    """
    if filetype == 'bigNarrowPeak':
        header = ('chr', 'start', 'end', 'name', 'score', 'strand',
                  'signalValue', 'pValue', 'qValue', 'peak')
        dtypes = {'end': int, 'start': int, 'score': float,
                  'signalValue': float, 'pValue': float, 'qValue': 'float',
                  'peak': float}
    else:
        raise NotImplementedError

    dtypes.update({k: str for k in header if k not in dtypes})

    return pd.read_table(filename, names=header, na_values='.', dtype=dtypes,
                         comment='#')


def filter_by_feature_type(feature, feature_type):
    if feature[2] == feature_type:
        return True


def add_chr(entry):
    entry.chrom = 'chr' + entry.chrom

    return entry


def count_reads_in_features(target, features_fn):
    """
    Callback function to count reads in features
    """
    return pybedtools.BedTool(target).intersect(
        b=features_fn,
        stream=True).count()


def get_sequences(bed, fasta):
    bt = pybedtools.BedTool(bed, from_string=True)
    bt.sequence(fi=fasta, name=True, s=True)

    with open(bt.seqfn) as f:
        seq = f.readlines()

    return [x.rstrip() for x in seq if not x.startswith('>')]


def plot_hist(values, log_x=False, title=""):
    t = "{2} distribution (μ={0:.2f}, σ={0:.2f})".format(
        np.mean(values), np.std(values), title)
    if log_x:
        p = figure(title=t, x_axis_type="log", )
        p.xaxis.axis_label = 'Log(x)'
    else:
        p = figure(title=t)
        p.xaxis.axis_label = 'x'

    hist, edges = np.histogram(values, density=True, bins='fd')

    p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
           fill_color="#036564", line_color="#036564")
    p.legend.location = "center_right"
    p.legend.background_fill_color = "darkgrey"
    p.yaxis.axis_label = 'Pr(x)'

    p.xaxis.major_label_text_font = 'helvetica'
    p.yaxis.major_label_text_font = 'helvetica'
    p.title.text_font = 'helvetica'
    p.xaxis.axis_label_text_font = 'helvetica'

    show(p)


def plot_vbar(values, title="counts", series=False):
    if not series:
        counts = pd.value_counts(values).to_frame(name='x_max')
    else:
        counts = values.to_frame(name='x_max')
    x_max = counts['x_max'].max()
    group = counts.sort_values(
        by='x_max', ascending=True).index.tolist()
    counts['x_min'] = 0
    counts = ColumnDataSource(counts.sort_values(by='x_max', ascending=False))

    p = figure(title=title, y_range=group, x_range=(0, x_max * 1.1))
    p.hbar(y="index", left='x_min', right='x_max', height=0.5, source=counts,
           fill_color="#036564", line_color="#036564")

    labels = LabelSet(x='x_max', y="index", text='x_max', level='glyph',
                      x_offset=5, y_offset=-5, source=counts,
                      render_mode='canvas',
                      text_font='helvetica', text_font_size='9pt')
    p.add_layout(labels)

    p.toolbar.active_drag = None

    p.ygrid.grid_line_color = None
    p.xaxis.axis_label = "Counts"
    p.yaxis.axis_line_color = None
    p.yaxis.major_tick_line_color = None
    p.outline_line_color = None

    p.xaxis.major_label_text_font = 'helvetica'
    p.yaxis.major_label_text_font = 'helvetica'
    p.title.text_font = 'helvetica'
    p.xaxis.axis_label_text_font = 'helvetica'
    show(p)


def count_by_feature(base_path):
    t_utr = pybedtools.BedTool(base_path + '3_utr.gff')
    f_utr = pybedtools.BedTool(base_path + '5_utr.gff')
    cds = pybedtools.BedTool(base_path + 'cds.gff')
    exon = pybedtools.BedTool(base_path + 'exon.gff')
    intergenic = pybedtools.BedTool(base_path + 'intergenic.bed')
    intron = pybedtools.BedTool(base_path + 'intron.gff')

    features = (t_utr, f_utr, cds, exon, intergenic, intron)
    pool = multiprocessing.Pool()
    results = pool.map(count_reads_in_features, features)
    counts = pd.Series(
        results, '3_utr 5_utr cds exon intergenic intron'.split())
    plot_vbar(counts, series=True, title='Peaks per feature')


def counts_by_ensembl_transcript_biotype(target, base_path):
    bed_obj = pybedtools.BedTool(target)
    ensb_gtf = pybedtools.BedTool(
        '/Volumes/biodb/genomes/homo_sapiens/GRCh38_90/GRCh38.90.gtf') \
        .filter(filter_by_feature_type, 'gene') \
        .each(add_chr) \
        .saveas()
    biotype_result = ensb_gtf.intersect(bed_obj, wa=True, wb=True)
    gene_type = {}
    for x in biotype_result:
        try:
            gene_type[x['gene_id']] = x['gene_biotype']
        except KeyError:
            pass
    plot_vbar(list(gene_type.values()), title='Peaks per gene biotype')
