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

$ git pull https://github.com/seppedewinter/gliftly.git
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

**Get aligned sequences across assemblies:**

```python

>>> from glifty.chain import ChainFile
>>> from glifty.tools import get_aligment, pretty_print_alignment
>>> cf = ChainFile.load_chain_file("hg19ToPanTro3.over.chain.gz")
>>> al = get_aligment(
        cf, 
        "chr2", 25_383_722, 25_383_722 + 1_000,
        "hg19.fa", "panTro3.fa"
    )
>>> pretty_print_alignment(*next(al))
source:	CTGTTATTTGACGGCTACGTATTTTTACTTTATTCACACAGTTTACATTC
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	CTGTTATTTGACGGCTACGTATTTTTACTTTATTCACACAGTTTACATTC
source:	AAAGTCAGAGGTGGATGTGAAATTTGAAAGGTTTTATTTCCTAACTACAG
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	AAAGTCAGAGGTGGATGTGAAATTTGAAAGGTTTTATTTCCTAACTACAG
source:	GCAGCTTTAAGAGGCTGATTATCTGCCACGACCCCCCAGGCTGGGAGGCG
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	GCAGCTTTAAGAGGCTGATTATCTGCCACGACCCCCCAGGCTGGGAGGCG
source:	GCAGCAGGGCAGGGGAGAGCAAGGGGCTTTGGGGTCGACCTCCTGGGGGA
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	GCAGCAGGGCAGGGGAGAGCAAGGGGCTTTGGGGTCGACCTCCTGGGGGA
source:	GGGTAGCCCTGGGGCCCCGCTGTGCCCTCACTCGCCCTTCTTGTAGGCGT
        |||||||||*|||||||*||||||||||||||||||||||||||||||||
target:	GGGTAGCCCCGGGGCCCGGCTGTGCCCTCACTCGCCCTTCTTGTAGGCGT
source:	TCTTGATGATGGCGTTTTTGAACAGCGTCACCAGGGGCGTCTGGCTCTTC
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	TCTTGATGATGGCGTTTTTGAACAGCGTCACCAGGGGCGTCTGGCTCTTC
source:	TCGGAGGTCATGAAACCGCCGTAGCGCTTGTCCTTGGGCGGGCTGCCCCA
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	TCGGAGGTCATGAAACCGCCGTAGCGCTTGTCCTTGGGCGGGCTGCCCCA
source:	GCGGAAGTGCTCCATCCTGTAGGGGCCCTCGTCCTTCTTCTCGGCCGCCA
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	GCGGAAGTGCTCCATCCTGTAGGGGCCCTCGTCCTTCTTCTCGGCCGCCA
source:	CCAGCAGGCTGTGCTCCAGGTCGGCCTGGGCCCCTGCGCCGTCATCGGCA
        ||||||||||||||||||||||||||||||||||*|||||||||||||||
target:	CCAGCAGGCTGTGCTCCAGGTCGGCCTGGGCCCCGGCGCCGTCATCGGCA
source:	GGGCCGTCGGGGCCATCTCCCTCCCGGAGTCGCTGGCCAGTCAGCTCCCT
        |||||||||||||||||||||||||||*||||||||||||||||||||||
target:	GGGCCGTCGGGGCCATCTCCCTCCCGGGGTCGCTGGCCAGTCAGCTCCCT
source:	CTTGAACTCCAGGGGGAAGGCCTCGGCCGACTCGTCCTCGGCGCCGTTAG
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	CTTGAACTCCAGGGGGAAGGCCTCGGCCGACTCGTCCTCGGCGCCGTTAG
source:	GGTACACCTTCACTGGGCGCCGCTTCTTGCCCACCGGCTTGCCCCAGCGG
        |||||||||||||*||||||||||||||||||||||||||||||||||||
target:	GGTACACCTTCACCGGGCGCCGCTTCTTGCCCACCGGCTTGCCCCAGCGG
source:	AAGTGCTCCATGGAGTAGGAGCGCTTGCCCTCGCGCGGGCCCGGCTTGGC
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	AAGTGCTCCATGGAGTAGGAGCGCTTGCCCTCGCGCGGGCCCGGCTTGGC
source:	ACCATCGCTGCGGGGCTCGGGGCCGCCCTCAGGCAGCGGGCCGCAGTCTT
        ||||||||||||||||||||||||||||||||||||||||||||*|||||
target:	ACCATCGCTGCGGGGCTCGGGGCCGCCCTCAGGCAGCGGGCCGCGGTCTT
source:	CGCCCGCTGAGACGTCCTCGCGCTTCTGCCCTGCgccgctgc---tgccg
        ||||||||||||||||||||||||||||||||||||||||||   |||*|
target:	CGCCCGCTGAGACGTCCTCGCGCTTCTGCCCTGCgccgctgccgctgctg
source:	ctgctgctgctgttgcGGCGGCCGAATCGGTCCCAGCGGAAGTGGCCCAT
        ||||||||||||||||||||||||||||||||||||||||||||||||||
target:	ctgctgctgctgttgcGGCGGCCGAATCGGTCCCAGCGGAAGTGGCCCAT
source:	GACGTACTTCCGGGGGTTCTCGGTCAGAGGCTGCTCGTCGCCATTTCCCG
        |||||||||||||||||||||||||||||||||||||||||||||*||||
target:	GACGTACTTCCGGGGGTTCTCGGTCAGAGGCTGCTCGTCGCCATTGCCCG
source:	GGAACATGGGAGTCTCGGCCGAGAGGTCGGGCTTGCAGGCCCGGATGCAC
        |||||||||||||||||||*||||||||||||||||||||||||||||||
target:	GGAACATGGGAGTCTCGGCGGAGAGGTCGGGCTTGCAGGCCCGGATGCAC
source:	TCCTGGGGGAAGACGCGAGGGCATGAGGGCAGCCCGTGCCCCGCACCCCG
        |||||||||||||||||||||||||||||||||||||||||*||||||||
target:	TCCTGGGGGAAGACGCGAGGGCATGAGGGCAGCCCGTGCCCTGCACCCCG
source:	GCCCGGCTGCCGCGCCCGTCACTGCGCCTAGGCCCTGGCCGCCCT-----
        |||||||||||||||||||||||||||||||||||||||||||||
target:	GCCCGGCTGCCGCGCCCGTCACTGCGCCTAGGCCCTGGCCGCCCTGGCCG
source:	----CGCCACGT
            ||||||||
target:	CCCTCGCCACGT

```

## Why use glifty over pyliftover?

Compared to [pyliftover](https://github.com/konstantint/pyliftover) Glifty makes use of the [NCLS](https://github.com/pyranges/ncls) datastructure for quick interval queries. Furthermore, whereas [pyliftover](https://github.com/konstantint/pyliftover) can only convert point intervals Glifty can convert intervals of arbitrary length.

