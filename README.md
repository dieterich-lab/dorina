[![Build Status](https://travis-ci.org/dieterich-lab/dorina.svg?branch=master)](https://travis-ci.org/dieterich-lab/dorina.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/dieterich-lab/dorina/badge.svg?branch=master)](https://coveralls.io/github/dieterich-lab/dorina?branch=master)

doRiNA
======

Database of posttranscriptional regulatory elements

Installation
------------

For doRiNA development, please use Python version > 3.4:

First, install bedtools, we recommend the latest release:
```
wget https://github.com/arq5x/bedtools2/archive/master.zip -O bedtools.zip
tar -xf bedtools.zip
cd bedtools2-master
make
```
Then, create a virtual environment and install doRiNA from source:
```
python3 -m venv dorina-dev
source activate dorina-dev
git clone https://github.com/dieterich-lab/dorina.git
cd dorina
pip install -r requirements.txt -r test_requirements.txt
pip install .
```

Quickstart
----------

TODO

Supported dataset
-----------------

## Regulatory
- RNA Binding Proteins obtained by cross-linking and
 immunoprecipitation derived experiments (CLIP)
- miRNA–RNA interactions obtained by clash experiments
- RNA binding predicted
- Single Nucleotide Variants stored in Ensembl


Data layout
-----------

At the moment, doRiNA 2 expects data to live in `/data/projects/doRiNA2`.

Within the data directory, doRiNA needs two directory structures, one for
genomes and one for regulators. Genomes are split into directories by clade,
species and assembly, e.g. `mammals/h_sapiens/hg19`. The assembly directory then
contains a number of GFF files with the different subsets of genomic data.
Creating these subsets from a GTF file downloaded from UCSC is accomplished by
running the `create_assembly.sh` script in the doRiNA source tree.

Regulators follow the identical split as with the genomes, just that they
additionally sort the regulators by miRNA and RBPs. Regulators are given as BED
files. To be picked up by doRiNA, they need to be accompanied by a JSON metadata
file with the same name but ending in .json instead of .bed.

As an example for the JSON format, take
`regulators/mammals/h_sapiens/hg19/RBP/PARCLIP_AGO1234_hg19.json`

```
{
    "summary": "AGO1-4 PAR-CLIP (Hafner 2010)*",
    "description": "Human Argonaute 1 through 4 (AGO1-4/EIF2C1-4) PAR-CLIP
clusters in HEK 293 cells",
    "methods": "AGO1_1, AGO2_1 and_2, AGO3_1 and_2, AGO4_1 and _2 PAR-CLIP cDNA
libraries (see supplement Table S2 in Hafner et al., 2010) were analyzed using
the PAR-CLIP Analysis Pipeline (PCP) (Jens et al. unpublished), essentially as
described in (Lebedeva et al. 2011). Clusters are regions of reference RNA
sequence contiguously covered by one or more reads from at least 2 different AGO
PAR-CLIP cDNA libraries.",
    "credits": "This track was produced by Marvin Jens in the N. Rajewsky lab at
the Berlin Institute of Medical Systems Biology at the Max Delbrück Center
Berlin. Marvin Jens, Zhou Fang, Sebastian Mackowiak, Andranik Ivanov and Jonas
Maaskola in the N. Rajewsky lab developed the computational pipeline for
analyzing PAR-CLIP experiments.",
    "references": [{
        "title": "Transcriptome-wide identification of RNA-binding protein and
microRNA target sites by PAR-CLIP.",
        "authors": [
            "Hafner M",
            "Landthaler M",
            "Burger L",
            "Khorshid M",
            "Hausser J",
            "Berninger P",
            "Rothballer A",
            "Ascano M Jr",
            "Jungkamp AC",
            "Munschauer M",
            "Ulrich A",
            "Wardle GS",
            "Dewell S",
            "Zavolan M",
            "Tuschl T"
        ],
        "pages": "141(1):129-41",
        "journal": "Cell",
        "year": "2010",
        "pubmed": "http://www.ncbi.nlm.nih.gov/pubmed/20371350"
    },
    {
        "title": "Transcriptome-wide Analysis of Regulatory Interactions of the
RNA-Binding Protein HuR.",
        "authors": [
            "Lebedeva S",
            "Jens M",
            "Theil K",
            "Schwanhäusser B",
            "Selbach M",
            "Landthaler M",
            "Rajewsky N"
        ],
        "pages": "43(3):340-52",
        "journal": "Mol. Cell",
        "year": "2011",
        "pubmed": "http://www.ncbi.nlm.nih.gov/pubmed/21723171"
    }]
}
```

License
-------

dorina is licensed under the GNU General Public Licence (GPL) version 3.
See `LICENSE` file for details.
