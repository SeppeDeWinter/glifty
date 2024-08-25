from glifty.chain import ChainFile
from typing import Iterator

def liftover(
        chain_file: ChainFile,
        chromosome: str,
        start: int,
        end: int,
        strand: str = "+"
    ) -> Iterator[tuple[str, int, int, str]]:
    """
    Perfrom UCSC style lifterover using a chain file

    Parameters
    ----------
    chain_file: ChainFile
        An instance of a ChainFile
    chromosome: str
        Query chromsome name
    start: int
        Query start position
    end: int
        Query end position
    strand: str
        Wether the query is on the positive ('+') or negative ('-') strand.

    Returns
    -------
    Iterator[tuple[str, int, int, str]]
        An iterator of lifted over coordinates in the form of 
        chromosome, start, end, strand.

    Example
    -------
    >>> from glifty.chain import ChainFile
    >>> from glifty.tools import liftover
    >>> cf = ChainFile.load_chain_file("hg19ToPanTro3.over.chain.gz")
    >>> lo = liftover(cf, "chr2", 25_383_722, 25_391_559)
    >>> next(lo)
    ('chr2A', 25494472, 25502382, '+')
    """
    results = chain_file.query(
        chrom = chromosome,
        start = start,
        end = end,
        strand = strand
    )
    for result in results:
        t_chrom, ((s_start, _, s_strand), (t_start, _, t_strand)) = result[0]
        _, ((_, s_end, _), (_, t_end, _)) = result[-1]
        if s_start != start or s_end != end or s_strand != strand:
            raise ValueError("Chainfile was not parsed correctly!")
        yield (t_chrom, t_start, t_end, t_strand)
