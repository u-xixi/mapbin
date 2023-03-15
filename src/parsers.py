import sys
from typing import Set, List, Tuple, Dict
from file_io import plain_reader
from collections import deque

def parse_pe_read_file_sep(r1_file, r2_file):
    pass


class NbRead:
    def __init__(self, header, seq):
        self.header = header
        self.seq = seq


def parse_contig_fa(contig_fasta:str) -> Dict[str, str]:
    """
    :param contig_fasta: contigs.fasta file path
    :return: {fasta header: seq}
    :return: {str: str}
    """
    contig_fa_lookup = {}

    def conclude_record(header, seq):
        contig_fa_lookup[header] = "".join(seq)

    with open(contig_fasta, 'r') as f:
        cur_header = ""
        cur_seq = []

        for iline in f:
            iline = iline.strip()
            if iline:
                if iline.startswith(">"):
                    if cur_header:
                        conclude_record(cur_header, cur_seq)
                    cur_header = iline.lstrip(">")
                    cur_seq = []
                else:
                    cur_seq.append(iline)
    return contig_fa_lookup


def parse_se_seqfile(file_content, mode: str) -> NbRead:
    # this returns items of an iterable

    ## for loop implementation: I don't know which is beter
    # if mode =='se_fa':
    #     cur_header, cur_seq = None, []
    #     for line in file_content:
    #         if line.startswith(">"):
    #             if cur_header:
    #                 yield NbRead(cur_header, cur_seq)
    #             cur_header, seq = line, []
    #         else:
    #             cur_seq.append(line)
    #     if cur_header:
    #         yield NbRead(cur_header, cur_seq)

    file_content = iter(file_content)
    if mode == "se_fa":
        cur_header, cur_seq = None, []
        while True:
            try:
                curline = next(file_content)
                if curline.startswith(">"):
                    if cur_header:
                        yield NbRead(cur_header, "".join(cur_seq))
                    cur_header, cur_seq = curline.strip()[1:], []
                else:
                    cur_seq.append(curline)
            except StopIteration:
                break
        if cur_header:
            yield NbRead(cur_header, cur_seq)


    elif mode == "se_fq":

        while True:
            try:
                curline = next(file_content)
                if curline.startswith("@"):
                    header = curline.strip()[1:]
                    seq = next(file_content).strip()
                    yield NbRead(header, seq)
                    next(file_content)
                    next(file_content)
            except StopIteration:
                break


def parse_pe_sep(r1_content: [], r2_content: [], mode: str):
    if len(r1_content) == len(r2_content):
        pair_content = zip(r1_content, r2_content)
        if mode == "pe_sep_fq":
            while True:
                try:
                    r1_line, r2_line = next(pair_content)
                    if r1_line.startswith("@") and \
                        r2_line.startswith("@"):
                        r1_header, r2_header = r1_line, r2_line
                        r1_seq, r2_seq = next(pair_content)
                        yield NbRead(r1_header, r1_seq), NbRead(r2_header, r2_seq)
                        next(pair_content)
                        next(pair_content)
                except StopIteration:
                    break


        elif mode == "pe_sep_fa":
            while True:
                try:
                    r1_line, r2_line = next(pair_content)
                    if r1_line.startswith(">") and \
                        r2_line.startswith(">"):
                        r1_header, r2_header = r1_line, r2_line
                        r1_seq, r2_seq = next(pair_content)
                        yield NbRead(r1_header, r1_seq), NbRead(r2_header, r2_seq)
                except StopIteration:
                    break


    else:
        print("[::Error::] invalid paired-end read files")
        # stub: return is None, and if that is the case, consider the seq files as non existent


def parse_pe_interleave(file_content: [], mode: str):
    if mode == "pe_il_fq":
        if len(file_content) % 8 == 0:
            file_content = iter(file_content)
            while True:
                try:
                    curline = next(file_content)
                    if curline.startswith("@"):
                        r1_header = curline.strip()[1:]
                        r1_seq = next(file_content)
                        next(file_content)
                        next(file_content)
                        next_seqrecord = next(file_content)
                        if next_seqrecord.startswith("@"):
                            r2_header = next_seqrecord.strip()[1:]
                            r2_seq = next(file_content).strip()
                        else:
                            sys.exit("[::Error::] invalid interleaved paired-end fastq.")
                except StopIteration:
                    break
        else:
            print("[::Warning::] ignoring invalid interleaved paired-end fastq.")

    elif mode == "pe_il_fa":
        if len(file_content)% 4 ==0:
            file_content = iter(file_content)
            while True:
                try:
                    curline = next(file_content)
                    if curline.startswith(">"):
                        r1_header = curline.strip()[1:]
                        r1_seq = next(file_content).strip()
                        next_seqrecord = next(file_content)
                        if next_seqrecord.startswith(">"):
                            r2_header = next_seqrecord.strip()[1:]
                            r2_seq = next(file_content).strip()
                            yield NbRead(r1_header, r1_seq), NbRead(r2_header, r2_seq)
                        else:
                            sys.exit("[::Error::] invalid interleaved paired-end fasta.")

                except StopIteration:
                    break
        else:
            print("[::Warning::] ignoring invalid interleaved paired-end fasta")



def parse_spades_contig_paths(all_contigs:Set[str], ctg_lookup: Dict, path_file:str) -> Dict[int, Set[str]]:
    """
    :param: ctg_lookup: Dict[str, NbContig]. contains contigs > min_ctg_len.
            keys are full fasta header
            values are formatted info to be used in infomap
    :param: contig_fa_lookup: Dict[str, str]. generated by network_utils.generate_ctg_lookup_alt()
    :returns: {contig: [(path_left_end, path_right_end)]}, contains only the contigs whose length > min_ctg_len
    :returns: {int: List[Tuple]}
    All overlap computing ignores the overlap lengths,
    because they are almost all the same, see gfa file's L lines"""

    def parse_path_line(path:[]) -> Set[str]:
        segs = set()
        for p in path:
            p_split = [i[:-1] for i in p.split(",")]
            segs.update(p_split)
            # ends.add(p_split[0])
            # ends.add(p_split[-1])
        return segs

    path_content = deque(plain_reader(path_file))
    path_segs = {}

    if path_content[0] in all_contigs \
            and path_content[1] not in all_contigs:
        # check if the format is correct by checking the first 2 lines
        is_rev = 0
        cur_header = path_content.popleft()  # initiate the loop with the first line
        cur_rev_header = cur_header + "\'"
        cur_path = []  # store path lines
        while path_content:
            iline = path_content.popleft()
            if iline in all_contigs:
                cur_header = iline
                is_rev = 0
                cur_rev_header = iline + "\'"

            elif iline == cur_rev_header:
                is_rev = 1
                # conclude the current entry
                if cur_path:
                    if cur_header in ctg_lookup:
                        cur_numid = ctg_lookup[cur_header].num_id
                        path_segs[cur_numid] = parse_path_line(cur_path)
                cur_path = []

            elif is_rev == 1:
                continue

            else:
                cur_path.append(iline.rstrip(";"))

    else:
        sys.exit("[::Error::] Incorrect paths file: {}".format(path_file))
    return path_segs


