from glifty.chain import ChainFile

from typing import Iterator
import pysam

GAP_CHAR = "-"

def liftover(
        chain_file: ChainFile,
        chromosome: str,
        start: int,
        end: int,
        strand: str = "+"
    ) -> Iterator[tuple[str, int, int, str]]:
    """
    Perfrom UCSC style lifterover using a chain file.

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

def _reverse_complement(seq: str) -> str:
    _complement = {
        ord("a"): ord("t"),
        ord("t"): ord("a"),
        ord("c"): ord("g"),
        ord("g"): ord("c"),
        ord("A"): ord("T"),
        ord("T"): ord("A"),
        ord("C"): ord("G"),
        ord("G"): ord("C"),
        ord("n"): ord("n"),
        ord("N"): ord("N")
    }
    return seq[::-1].translate(_complement)

def _get_sequence(
        genome: pysam.FastaFile,
        chromosome: str,
        start: int,
        end: int,
        strand: str
) -> str:
    if strand == "+":
        seq = str(genome[chromosome][start:end])
    elif strand == "-":
        chrom_size = genome.get_reference_length(chromosome)
        seq = str(
            genome[chromosome][
                chrom_size - end: chrom_size - start])
        seq = _reverse_complement(seq)
    else:
        raise ValueError(f"Invalid strand: {strand}")
    return seq

def get_alignment(
        chain_file: ChainFile,
        chromosome: str,
        start: int,
        end: int,
        genome_filename_source: str,
        genome_filename_target: str,
        strand: str = "+"
) -> Iterator[tuple[str, str]]:
    """
    Get aligned sequences from chain file.

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
    genome_filename_source: str
        filename of the source genome fasta
    genome_filename_target: str
        filename of the target genome fasta
    strand: str
        Wether the query is on the positive ('+') or negative ('-') strand.

    Returns
    -------
    Iterator[tuple[str, str]]
        An iterator of aligned sequences

    Example
    -------
    >>> from glifty.chain import ChainFile
    >>> from glifty.tools import get_aligment
    >>> cf = ChainFile.load_chain_file("hg19ToPanTro3.over.chain.gz")
    >>> al = get_aligment(
            cf, 
            "chr2", 25_383_722, 25_391_559,
            hg19.fa, panTro3.fa
        )
    >>> # TODO
    """
    results = chain_file.query(
        chrom = chromosome,
        start = start,
        end = end,
        strand = strand
    )
    source_genome = pysam.FastaFile(genome_filename_source)
    target_genome = pysam.FastaFile(genome_filename_target)
    for result in results:
        # load the whole DNA sequence without gaps.
        t_chrom, ((s_start_i, _, s_strand), (t_start_i, _, t_strand)) = result[0]
        _, ((_, s_end_i, _), (_, t_end_i, _)) = result[-1]
        source_seq = _get_sequence(
            genome = source_genome,
            chromosome = chromosome,
            start = s_start_i,
            end = s_end_i,
            strand = s_strand
        )
        target_seq = _get_sequence(
            genome = target_genome,
            chromosome = t_chrom,
            start = t_start_i,
            end = t_end_i,
            strand = t_strand
        )
        # offset used to get relative indices in the loaded DNA sequences
        s_offset = s_start_i
        t_offset = t_start_i
        # s_end_prev and t_end_prev stores the end coordinates of the previous
        # iteration of the for loop, this to calculate the number of gaps in the
        # alignment.
        _, ((s_start, s_end_prev, _), (t_start, t_end_prev, _)) = result[0]
        source_seq_aligned = source_seq[
            s_start - s_offset: s_end_prev - s_offset
        ]
        target_seq_aligned = target_seq[
            t_start - t_offset: t_end_prev - t_offset
        ]
        for _, ((s_start, s_end, _), (t_start, t_end, _)) in result[1:]:
            # Add gaps
            source_seq_aligned += source_seq[
                s_end_prev - s_offset: s_start - s_offset
            ]
            target_seq_aligned += (s_start - s_end_prev) * GAP_CHAR

            target_seq_aligned += target_seq[
                t_end_prev - t_offset: t_start - t_offset
            ]
            source_seq_aligned += (t_start - t_end_prev) * GAP_CHAR

            # Add aligned block
            source_seq_aligned += source_seq[
                s_start - s_offset: s_end - s_offset
            ]
            target_seq_aligned += target_seq[
                t_start - t_offset: t_end - t_offset
            ]
            s_end_prev = s_end
            t_end_prev = t_end

        yield (source_seq_aligned, target_seq_aligned)

    source_genome.close()
    target_genome.close()
