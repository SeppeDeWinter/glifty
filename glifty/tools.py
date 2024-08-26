from glifty.chain import ChainFile

from typing import Iterator, Sequence, Any, TypeVar
import pysam

SEQ_OR_STR = TypeVar("SEQ_OR_STR", str, Sequence[Any])
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

def pretty_print_alignment(source: str, target: str, width = 50):
    for i in range(0, max(len(source), len(target)), width):
        print(f"source:\t{source[i:i+width]}")
        print( "       \t", end = "")
        for s, t in zip(source[i:i+width], target[i:i+width]):
            if s == GAP_CHAR or t == GAP_CHAR:
                print(" ", end = "")
            elif s == t:
                print("|", end = "")
            else:
                print("*", end = "")
        print()
        print(f"target:\t{target[i:i+width]}")
        print()

def _align_sequence(
        source_seq: SEQ_OR_STR,
        target_seq: SEQ_OR_STR,
        lo_result: list[tuple[str, tuple[tuple[int, int, str], tuple[int, int, str]]]],
        GAP: SEQ_OR_STR = GAP_CHAR
)  -> tuple[list[SEQ_OR_STR], list[SEQ_OR_STR]]:
    source_aligned = []
    target_aligned = []
    _, ((s_offset, _, _), (t_offset, _, _)) = lo_result[0]
    # s_end_prev and t_end_prev stores the end coordinates of the previous
    # iteration of the for loop, this to calculate the number of gaps in the
    # alignment.
    _, ((s_start, s_end_prev, _), (t_start, t_end_prev, _)) = lo_result[0]
    source_aligned.append(source_seq[
        s_start - s_offset: s_end_prev - s_offset
    ])
    target_aligned.append(target_seq[
        t_start - t_offset: t_end_prev - t_offset
    ])
    for _, ((s_start, s_end, _), (t_start, t_end, _)) in lo_result[1:]:
        # Add gaps
        source_aligned.append(source_seq[
            s_end_prev - s_offset: s_start - s_offset
        ])
        for _ in range(s_start - s_end_prev):
            target_aligned.append(GAP)

        target_aligned.append(target_seq[
            t_end_prev - t_offset: t_start - t_offset
        ])
        for _ in range(t_start - t_end_prev):
            source_aligned.append(GAP)

        # Add aligned block
        source_aligned.append(source_seq[
            s_start - s_offset: s_end - s_offset
        ])
        target_aligned.append(target_seq[
            t_start - t_offset: t_end - t_offset
        ])
        s_end_prev = s_end
        t_end_prev = t_end
    return source_aligned, target_aligned

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
        source_seq_aligned_, target_seq_aligned_ = _align_sequence(
            source_seq, target_seq,
            result, GAP_CHAR
        )

        source_seq_aligned = "".join(source_seq_aligned_)
        target_seq_aligned = "".join(target_seq_aligned_)

        assert source_seq_aligned.replace(GAP_CHAR, "") == source_seq, "Aligned sequence does not macth source sequence..."

        assert target_seq_aligned.replace(GAP_CHAR, "") == target_seq, "Aligned sequence does not macth target sequence..."

        yield (source_seq_aligned, target_seq_aligned)

    source_genome.close()
    target_genome.close()
