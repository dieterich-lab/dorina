[![Build Status](https://travis-ci.org/dieterich-lab/dorina.svg?branch=master)](https://travis-ci.org/dieterich-lab/dorina.svg?branch=master)
[![Coverage Status](https://coveralls.io/repos/github/dieterich-lab/dorina/badge.svg?branch=master)](https://coveralls.io/github/dieterich-lab/dorina?branch=master)

doRiNA
======

Database of posttranscriptional regulatory elements and function

Installation
------------

First, install bedtools, we recommend the latest release:
```bash
wget https://github.com/arq5x/bedtools2/archive/master.zip -O bedtools.zip
tar -xf bedtools.zip
cd bedtools2-master
make
```
Then, create a virtual environment and install doRiNA from source:
```bash
python3 -m venv dorina-dev
source activate dorina-dev
git clone https://github.com/dieterich-lab/dorina.git
cd dorina
pip install -r requirements.txt -r test_requirements.txt
pip install .
```

For doRiNA development, please use Python version > 3.4.

Usage
------

Dorina requires the files to be setup locally. Let `/path/to/datasets/` be the the setup path. 
`/path/to/datasets/` should contain :
- `/path/to/datasets/regulators`
- `/path/to/datasets/genomes`

Both sub-directories have the same structure:
- `/path/to/datasets/genomes/{organism}/{assembly}/`
- `/path/to/datasets/regulators/{organism}/{assembly}/`

Regulators can be obtained from (here)[http://dorina.mdc-berlin.de/regulators]. Each bed file should be accompanied with metadata file with the same name, but the `.json` extension:

```json
{ "id": "d_melanogaster", "label": "Fly", "scientific": "Drosophila melanogaster", "weight": 3}
```

The genome directory contains genome annotation separated into several regions for further filtering.

Given the directory structure is correct, the following command should retrieve the miR-1247 regulators of human genome assembly `hg19`:

```bash
dorina run 'hg19' --seta 'hsa-miR-1247|CLASH' -p /path/to/datasets/ > miR-1247.bed
```

To list the avaiable data sources, use:
```bash
dorina genomes -p /path/to/datasets/ | less 
dorina regulators -p /path/to/datasets/ | less 
```

Supported dataset
-----------------

## Regulatory
- RNA Binding Proteins obtained by cross-linking and
 immunoprecipitation derived experiments (CLIP)
- miRNA–RNA interactions obtained by clash experiments
- Predicted RNA interactions

## Expression
Dorina can filter the regulatory datasets by expression

## Variation
Retrieves single nucleotide variants co-occurring with regulatory elements 

License
-------

dorina is licensed under the GNU General Public Licence (GPL) version 3.
See `LICENSE` file for details.
