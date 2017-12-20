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
- Predicted RNA interactions

Data layout
-----------


License
-------

dorina is licensed under the GNU General Public Licence (GPL) version 3.
See `LICENSE` file for details.
