# GLIFTY: Genomic LIFTtover in pYthon

A Python library for quick conversion of genomic coordinates between assemblies.

## Install

**PyPi**

```bash

$ pip install glifty

```

**From Source**

Using [poetry](https://python-poetry.org/)

```bash

$ git clone https://github.com/SeppeDeWinter/glifty.git
$ cd glifty
$ poetry install

```

## Usage

**Liftover of coordinates between assemblies:**

```python

>>> from glifty.chain import ChainFile
>>> from glifty.tools import liftover
>>> cf = ChainFile.load_chain_file("hg19ToPanTro3.over.chain.gz")
>>> lo = liftover(cf, "chr2", 25_383_722, 25_391_559)
>>> next(lo)
('chr2A', 25494472, 25502382, '+')

```

## Why use glifty over pyliftover?

Compared to [pyliftover](https://github.com/konstantint/pyliftover) Glifty makes use of the [NCLS](https://github.com/pyranges/ncls) datastructure for quick interval queries. Furthermore, whereas [pyliftover](https://github.com/konstantint/pyliftover) can only convert point intervals Glifty can convert intervals of arbitrary length.

