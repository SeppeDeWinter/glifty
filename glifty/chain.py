from ncls import NCLS
import gzip
import io
from dataclasses import dataclass
from typing import Iterator
from typing_extensions import Self
from tqdm import tqdm

@dataclass
class ChainHeader:
    score: int
    sName: str
    sSize: int
    sStrand: str
    sStart: int
    sEnd: int
    tName: str
    tSize: int
    tStrand: str
    tStart: int
    tEnd: int
    _id: int

    @classmethod
    def from_header_string(cls, header: str) -> Self:
        try:
            (
                chain,
                score,
                sName,
                sSize,
                sStrand,
                sStart,
                sEnd,
                tName,
                tSize,
                tStrand,
                tStart,
                tEnd,
                _id 
            ) = header.strip().split()
            
            if chain != "chain":
                raise ValueError(f"Invalid header: {header}")
            if sStrand not in ["+", "-"] or tStrand not in ["+", "-"]:
                raise ValueError(f"Invalid header: {header}")

            return cls(
                score = int(score),
                sName = sName,
                sSize = int(sSize),
                sStrand = sStrand,
                sStart = int(sStart),
                sEnd = int(sEnd),
                tName = tName,
                tSize = int(tSize),
                tStrand = tStrand,
                tStart = int(tStart),
                tEnd = int(tEnd),
                _id  = int(_id)
            )
        except ValueError:
            raise ValueError(f"Invalid header: {header}")

class Chain:
    def _readline(self, file: gzip.GzipFile | io.TextIOWrapper) -> str:
        data = file.readline()
        if isinstance(data, bytes):
            return data.decode().strip()
        else:
            return data.strip()

    def __init__(
        self,
        header: ChainHeader,
        file: gzip.GzipFile | io.TextIOWrapper,
        show_progress: bool = False
    ):
        self.header = header
        fields = self._readline(file).split()
        # lists to store start/end coordinates of the source (target) genome
        s_starts: list[int] = []
        s_ends: list[int] = []
        # lists to store start/end coordinates of the target (query) genome
        t_starts: list[int] = []
        t_ends: list[int] = []
        # start of the alignment block for source and target genome
        s_from = header.sStart
        t_from = header.tStart
        if show_progress:
            pbar = tqdm(
                total = header.sEnd - header.sStart,
                unit = "bps",
                leave = False,
                desc = header.sName)
        while len(fields) == 3:
            try:
                # size: number of aligned bases without gaps
                # dt: number of gaps in the source genome (after the aligned section)
                # dq: number of gaps in the target genome (after the aligned section)
                size, dt, dq = map(int, fields)
                s_starts.append(s_from)
                s_ends.append(s_from + size)
                t_starts.append(t_from)
                t_ends.append(t_from + size)
                s_from += size + dt
                t_from += size + dq
                if show_progress:
                    pbar.update(size + dt)
                fields = self._readline(file).split()
            except ValueError:
                raise ValueError(f"Invalid field:{' '.join(fields)}")
        if len(fields) != 1:
            raise ValueError(
                f"Expected a single number on the last line, found {' '.join(fields)}.")
        try:
            size = int(fields[0])
            s_starts.append(s_from)
            s_ends.append(s_from + size)
            t_starts.append(t_from)
            t_ends.append(t_from + size)
        except ValueError:
            raise ValueError(f"Invalid field:{' '.join(fields)}")

        if s_ends[-1] != header.sEnd or t_ends[-1] != header.tEnd:
            raise ValueError(
                f"chain data is inconsistent with alignment size specified in header!\n{header}")

        # Create Nested containment list for easy querying.
        # ids will be used to retrieve query coordinates later
        self._ncls = NCLS(
            starts = s_starts,
            ends = s_ends,
            ids = list(range(len(s_starts)))
        )
        self._t_coords = list(zip(t_starts, t_ends))

    def query(
            self, start: int, end: int
    ) -> Iterator[tuple[str, tuple[tuple[int, int, str], tuple[int, int, str]]]]:
        overlaps = list(self._ncls.find_overlap(start, end))
        for i, (s_start, s_end, idx) in enumerate(overlaps):
            t_start, t_end = self._t_coords[idx]
            # start and end positions need to be adjusted
            # the start in the chain file might be upstream of the query start
            # similarly the end might be downstrean of the query end
            if i == 0:
                if s_start < start:
                    t_start = t_start + (start - s_start)
                    s_start = start

            if i == len(overlaps) - 1:
                if s_end > end:
                    t_end = t_end - (s_end - end)
                    s_end = end

            yield (
                self.header.tName, 
                (
                    (s_start, s_end, self.header.sStrand),
                    (t_start, t_end, self.header.tStrand)
                )
            )

@dataclass
class ChainFile:
    chains: dict[tuple[str, str], list[Chain]]
    filename: str

    def __repr__(self) -> str:
        return f"{self.filename} containing {len(self.chains)} chains."


    @classmethod
    def load_chain_file(cls, filename: str):
        if filename.endswith(".gz"):
            file: gzip.GzipFile | io.TextIOWrapper = gzip.open(filename)
        else:
            file = open(filename)
        chains: dict[tuple[str, str], list[Chain]] = {}
        header_line = file.readline()
        while header_line:
            header = ChainHeader.from_header_string(
                header_line.decode() if isinstance(header_line, bytes) else header_line
            )
            if (header.sName, header.sStrand) in chains.keys():
                chains[(header.sName, header.sStrand)].append(Chain(header, file))
            else:
                chains[(header.sName, header.sStrand)] = [Chain(header, file)]
            _ = file.readline() # there is an empty line after each chain
            header_line = file.readline()
        return cls(chains = chains, filename = filename)

    def query(
        self, chrom:str, start: int, end: int, strand: str = "+"
    ) -> list[list[tuple[str, tuple[tuple[int, int, str], tuple[int, int, str]]]]]:
        if (chrom, strand) not in self.chains.keys():
            raise ValueError(f"{(chrom, strand)} not in chain file!")
        overlaps = []
        for c in self.chains[(chrom, strand)]:
            results = list(c.query(start, end))
            if len(results) > 0:
                overlaps.append(results)
        return overlaps
