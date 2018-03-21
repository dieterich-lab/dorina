#!/usr/bin/env python
# -*- coding: utf-8
"""
Created on 16:54 19/03/2018 2018 

"""
import pandas as pd
from bokeh.io import curdoc
from bokeh.layouts import layout, widgetbox
from bokeh.models import ColumnDataSource, HoverTool, Div
from bokeh.models.widgets import Slider, Select, TextInput
from bokeh.plotting import figure

from os.path import dirname, join

deseq = pd.read_csv(open(join(dirname(__file__), 'test_deseq_out.csv')))
deseq = deseq.loc[:10000]
deseq = deseq.rename(columns={'Unnamed: 0': 'Gene name'})

axis_map = {
    "Mean of normalized counts": "baseMean",
    "Fold change": "log2FoldChange",
    "Adjusted pvalue": 'padj',
    'p-value': 'pvalue'
}


def select_genes():
    selected = deseq.copy()
    gene_val = gene.value.strip()
    selected = selected[
        (selected['padj'] < padj.value)
    ]
    if gene_val:
        selected = selected[
            selected['Gene name'].str.contains(gene_val) == True]
    return selected


def update():
    df = select_genes()

    x_name = axis_map[x_axis.value]
    y_name = axis_map[y_axis.value]
    p.title.text = "%d genes selected" % len(df)
    p.xaxis.axis_label = x_axis.value
    p.yaxis.axis_label = y_axis.value
    source.data = dict(
        x=df[x_name],
        y=df[y_name],
        name=df["Gene name"],
    )


# Input controls
x_axis = Select(title="X Axis", options=sorted(axis_map.keys()),
                value="Adjusted pvalue")
y_axis = Select(title="Y Axis", options=sorted(axis_map.keys()),
                value="Fold change")
gene = TextInput(title="Gene name contains")
padj = Slider(title="Adjusted pval", value=0.05, start=deseq['padj'].min(),
              end=deseq['padj'].max(), step=0.05)

source = ColumnDataSource(data=dict(x=[], y=[], name=[]))

hover = HoverTool(tooltips=[
    ("Counts", "@x"),
    ("FC", "@y"),
    ("Name", "@name")
])

p = figure(plot_height=600, plot_width=700, title="", toolbar_location=None,
           tools=[hover])
p.circle(x="x", y="y", source=source, size=7, line_color=None)

controls = [gene, x_axis, y_axis, padj]
for control in controls:
    control.on_change('value', lambda attr, old, new: update())

sizing_mode = 'fixed'  # 'scale_width' also looks nice with this example

inputs = widgetbox(*controls, sizing_mode=sizing_mode)

desc = Div(text=open(join(dirname(__file__), "description.html")).read(),
           width=800)

l = layout([
    [desc],
    [inputs, p],
], sizing_mode=sizing_mode)

update()  # initial load of the data

curdoc().add_root(l)
