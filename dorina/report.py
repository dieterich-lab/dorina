#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 13:48 16/04/2018 2018 

"""
import multiprocessing
import pathlib

import numpy as np
import pandas as pd
import pybedtools
from bokeh.layouts import gridplot
from bokeh.models import ColumnDataSource, LabelSet
from bokeh.plotting import figure, output_file


def featuretype_filter(feature, featuretype):
    if feature[2] == featuretype:
        return True


def count_reads_in_features(bed, features_fn):
    """
    Callback function to count reads in features
    """
    return bed.intersect(b=features_fn, stream=True).count()


def get_sequences(bed, fasta):
    bed.sequence(fi=fasta, name=True, s=True)
    with open(bed.seqfn) as f:
        seq = f.readlines()

    return [x.rstrip() for x in seq if not x.startswith('>')]


def plot_hist(values, logx=False, title=""):
    t = "{2} distribution (μ={0:.2f}, σ={0:.2f})".format(
        np.mean(values), np.std(values), title)
    if logx:
        p1 = figure(title=t, x_axis_type="log")
        p1.xaxis.axis_label = 'Log(x)'
    else:
        p1 = figure(title=t)
        p1.xaxis.axis_label = 'x'

    hist, edges = np.histogram(values, density=True, bins='fd')

    p1.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:],
            fill_color="#036564", line_color="#036564")
    p1.legend.location = "center_right"
    p1.legend.background_fill_color = "darkgrey"
    p1.yaxis.axis_label = 'Pr(x)'

    p1.xaxis.major_label_text_font = 'helvetica'
    p1.yaxis.major_label_text_font = 'helvetica'
    p1.title.text_font = 'helvetica'
    p1.xaxis.axis_label_text_font = 'helvetica'

    return p1


def plot_vbar(values, title="counts", series=False, count=False, keys=None):
    if count:
        counts = pd.value_counts(values).to_frame(name='x_max')
    else:
        counts = values.to_frame(name='x_max')

    if keys:
        counts = counts.loc[keys]
    else:
        counts = counts.sort_values(by='x_max', ascending=True)

    x_max = counts['x_max'].max()
    group = counts.index.tolist()
    counts['x_min'] = 0
    counts = ColumnDataSource(counts)
    #     hover = HoverTool(tooltips=[( 'Counts', '@x_max')])
    p = figure(title=title, y_range=group, x_range=(0, x_max * 1.1),
               plot_width=500, plot_height=750)
    p.hbar(y="index", left='x_min', right='x_max', height=0.5,
           source=counts, fill_color="#036564", line_color="#036564")
    labels = LabelSet(x='x_max', y="index", text='x_max',
                      level='glyph', x_offset=5, y_offset=-7.5,
                      source=counts, render_mode='canvas',
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

    return p


def add_chr(entry):
    entry.chrom = 'chr' + entry.chrom

    return entry


def plot_chr_counts(assembly, dataframe):
    chromsizes = {k: pybedtools.chromsizes(assembly)[k][1] - \
                     pybedtools.chromsizes(assembly)[k][0]
                  for k in pybedtools.chromsizes(assembly)}

    keys = dataframe['chrom'].value_counts(
    ).sort_values(ascending=True).index.tolist()
    return gridplot([[
        plot_vbar(pd.Series(dataframe['chrom']), count=True, keys=keys,
                  title='Counts per chromossome'),
        plot_vbar(pd.Series(chromsizes), keys=keys,
                  title=assembly + ' Chromossome size')]])


def plot_feat_counts(bt, datadir, n_proc=1):
    def count_reads_in_features_this(features):
        return count_reads_in_features(bt, features_fn=features)

    def total_feature_length(bed_Obj):
        df = bed_Obj.to_dataframe()
        return sum(df['end'] - df['start'])

    t_utr = pybedtools.BedTool(datadir + '/3_utr.gff')
    f_utr = pybedtools.BedTool(datadir + '/5_utr.gff')
    cds = pybedtools.BedTool(datadir + '/cds.gff')
    exon = pybedtools.BedTool(datadir + '/exon.gff')
    intergenic = pybedtools.BedTool(datadir + '/intergenic.bed')
    intron = pybedtools.BedTool(datadir + '/intron.gff')

    features = (t_utr, f_utr, cds, exon, intergenic, intron)
    feat_names = '3_utr 5_utr cds exon intergenic intron'.split()

    with multiprocessing.Pool(processes=n_proc) as pool:
        results = pool.map(count_reads_in_features_this, features)

    with multiprocessing.Pool(processes=n_proc) as pool:
        features_length = pool.map(total_feature_length, features)

    counts_per_feature = pd.Series(results, feat_names)

    features_length = pd.Series(features_length, feat_names)

    return gridplot([[
        plot_vbar(counts_per_feature, series=True, title='Peaks per feature'),
        plot_vbar(features_length, series=True, title='Feature length')]])


def plot_biotype_counts(bt, ensembl_gtf):
    bt_gtf = pybedtools.BedTool(
        ensembl_gtf) \
        .filter(featuretype_filter, 'gene') \
        .each(add_chr) \
        .saveas()
    biotype_result = bt_gtf.intersect(bt, wa=True, wb=True)

    gene_type = {}
    for x in biotype_result:
        try:
            gene_type[x['gene_id']] = x['gene_biotype']
        except KeyError:
            pass

    return plot_vbar(list(gene_type.values()), title='Peaks per gene biotype')


def main(target, regulator=None, fasta=None, output_dir=None,
         assembly='hg38', datadir=None, n_proc=1, ensembl_gtf=None):
    if output_dir is None:
        output_dir = pathlib.Path.cwd()

    bt = pybedtools.BedTool(target)
    if regulator:
        bt = bt.filter(lambda x: regulator in x.name).saveas()

    df = bt.to_dataframe()
    if fasta:
        df['seq'] = get_sequences(bt, fasta=fasta)

    output_file(
        str(output_dir / 'score_dist.png'),
        plot_hist(df['score'], logx=True, title='Score'))

    output_file(
        str(output_dir / 'peak_length.png'),
        plot_hist(df['end'] - df['start'], title='Peak length'))

    output_file(
        str(output_dir / 'count_per_chr.png'),
        plot_chr_counts(assembly, df))

    output_file(
        str(output_dir / 'count_per_feature.png'),
        plot_feat_counts(bt, datadir, n_proc=n_proc))

    output_file(
        str(output_dir / 'count_per_biotype.png'),
        plot_biotype_counts(bt, ensembl_gtf))

    return 0