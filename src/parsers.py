from typing import Set, List, Tuple, Dict, Iterable
import glob
from pysam import AlignmentFile,AlignedSegment
from seq_objects import *
import gzip
import time
import multiprocessing as mp
import os
from operator import add
from collections import Counter, namedtuple

from random import sample


def isgzip_reader(file_path):
    try:
        if file_path.endswith(".gz") or \
                file_path.endswith("gzip"):
            return gzip.open(file_path, "rt")

        else:
            return open(file_path, 'r')

    except ValueError:
        print("[=== Error ===] File corrupted or file path invalid!")


def parse_contig_fa(contig_fasta:str) -> Dict[str, str]:
    """
    :param contig_fasta: contigs.fasta file path
    :return: {fasta header: seq}
    :return: {str: str}
    """
    contig_fa_lookup = {}

    def conclude_record(header, seq):
        contig_fa_lookup[header] = "".join(seq)

    f = isgzip_reader(contig_fasta)
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
    conclude_record(cur_header, cur_seq)
    f.close()
    return contig_fa_lookup


def make_binned_lookup(args, ctg_lookup: Dict[str, NbContig]):
    # def make_binned_lookup(bin_dir: str, bin_suffix: str, ctg_lookup: Dict[str, NbContig]):
    bin_lookup = {}  # ctg_id: bin_idx
    i = 1
    binned_files = glob.glob(args.bin_dir + "*." + args.bin_suffix)

    if binned_files:
        for bin_fa in binned_files:
            print("[=== Info ===] Read contig fasta: {}".format(bin_fa))
            with open(bin_fa, 'r') as f:
                for iline in f:
                    if iline.startswith(">"):
                        ctg_name = iline.rstrip()[1:].split(" ")[0]
                        if ctg_name in ctg_lookup:
                            bin_lookup[ctg_name] = i
            i += 1

        with open(args.init_clu, 'w') as fw:
            fw.write("# node module flow\n")
            for k, v in bin_lookup.items():
                fw.write("{} {}\n".format(ctg_lookup[k].num_id, v))

    else:
        print("[=== Warning ===] Find no sequence files in "
              "directory \"{}\". \"Binned\" mode off".format(args.bin_dir + "*." + args.bin_suffix))

    return bin_lookup


class LinkedCtg():
    # contig is the basic unit of this network
    # use this class as the basis
    # it might help reduce the overhead
    # and allows modifying or adding info to the contig
    __slots__ = ['id', 'interval_len', 'cov', 'binid', 'direction']

    def __init__(self, id:str, interval_len=NotImplemented,
                 cov=NotImplemented, binid=NotImplemented,
                 direction=NotImplemented):
        """
        :param id: NODE_1
        :param nbases: linked bp
        :param cov: the coverage of the linked part
        """
        self.id = id
        self.interval_len = interval_len
        self.cov = cov
        self.binid = binid
        self.direction = direction

class ReadAlignment():
    # contig is the basic unit of this network
    # use this class as the basis
    # it might help reduce the overhead
    # and allows modifying or adding info to the contig
    __slots__ = ['t_start', 't_end']

    def __init__(self, t_start:int, t_end:int):
        """
        :param id: NODE_1
        :param nbases: linked bp
        :param cov: the coverage of the linked part
        """
        self.t_start = t_start
        self.t_end = t_end


def get_n_subprocesses(n_threads: int, n_files:int):
    n_sys = mp.cpu_count()
    nproc = min(n_sys, n_threads, n_files)
    if n_sys < n_threads:
        print("[=== Warning ===] Detected only {} CPUs instead of {}. "
              "Set number of CPUs to {}".format(n_sys, n_threads, nproc))
    return nproc


# N_PROCESS = get_n_subprocesses


def process_bams(aln_files: List[str], aln_query_perc: float,
                 ctg_lookup: Dict[str, NbContig],
                 nprocess: int,
                 bin_lookup: Dict[str, int],
                 rlen: int,
                 max_insert_size: int):

    def update_ctg_links(single_links: Dict[Tuple[str, str], int]):
        for k, v in single_links.items():
            if k in ctg_links:
                ctg_links[k] += v
            else:
                ctg_links[k] = v

    # check files:
    aln_files = set(aln_files)
    print("[=== Info ===] {} unique bam files are given".format(len(aln_files)))
    for af in aln_files:
        if not os.path.isfile(af):
            raise FileNotFoundError(af)

    ctg_links = {}

    reports = []
    with mp.Pool(processes=nprocess) as pool:
        for af in aln_files:
            args_single = (af, aln_query_perc, ctg_lookup, bin_lookup, rlen, max_insert_size)
            reports.append(pool.apply_async(process_a_bam,
                                            args_single, ))

        success_all = False
        fail_some = False

        while not (success_all or fail_some):
            time.sleep(5)
            success_all = all(tr.ready() and tr.successful() for tr in reports)
            fail_some = any(tr.ready() and not tr.successful() for tr in reports)

            if success_all:
                pool.close()
            if fail_some:
                pool.terminate()

        pool.join()

    # exit the program and inform the user
    for af, r in zip(aln_files, reports):
        if r.ready() and not r.successful():
            print("[=== Error ===] An error occurred during the parallel processing of {}. "
                  "Program terminated. Please check the following error message".format(af))
            print(r.get())

    # if all successful, proceed:
    for r in reports:
        bam_ctg_links = r.get()
        if not ctg_links:
            print("======= parser line 273, not all ctg links, len(bam ctg link)", len(bam_ctg_links))
            ctg_links = bam_ctg_links
        else:
            print("======= parser line 276, all_ctg_links, len(bam ctg link)", len(bam_ctg_links))
            update_ctg_links(bam_ctg_links)

    return ctg_links


def process_a_bam(aln_file, aln_query_perc, ctg_lookup,
                  bin_lookup: Dict[str, int],
                  rlen: int,
                  max_insert_size: int
                  ) -> Dict[Tuple[str, str], int]:

    def get_cov_intervals(read_pairs: List[ReadAlignment], max_cov_gap: float):
        # get the intervals where no read cover over
        read_pairs = sorted(read_pairs, key=lambda x: x.t_start)
        print("==========read locs:", [(i.t_start, i.t_end) for i in read_pairs])
        prev_read = read_pairs[0]
        valid_intervals = [prev_read.t_start, prev_read.t_end]  # start, end
        covered_bases = valid_intervals[1] - valid_intervals[0]

        for i in range(1, len(read_pairs)):
            cur_read = read_pairs[i]

            if cur_read.t_start < valid_intervals[1] + max_cov_gap:
                covered_bases += cur_read.t_end - cur_read.t_start
                valid_intervals[1] = cur_read.t_end
            else:
                return [], 0

        return valid_intervals, covered_bases

    def _check_end_poz(start_poz: int, end_poz: int, template_len: int):
        # if the distance to the ends is above 1000, call it off
        till_template_end_len = template_len - end_poz
        if start_poz < till_template_end_len:
            return start_poz, 1

        else:
            return till_template_end_len, -1

    def _get_bin_id(ctg_id: str):
        if ctg_id in bin_lookup:
            return bin_lookup[ctg_id]
        else:
            return -1

    def _get_interval_len(_interval):
        return _interval[1] - _interval[0]

    def _calc_cov(_interval_len, _nbases:int):
        return _nbases/_interval_len

    def _update_ctg2ctgs(_ctg2ctgs, _ctg, _il_ctg, _interval_len, _cov, _binid, _direction):
        _linkedctg = LinkedCtg(_il_ctg, _interval_len, _cov, _binid, _direction)
        if _ctg in _ctg2ctgs:
            _ctg2ctgs[_ctg].append(_linkedctg)

        else:
            _ctg2ctgs[_ctg] = [_linkedctg]

    def filter_links(_raw_links: Dict[Tuple[str, str], List[Tuple[ReadAlignment, ReadAlignment]]]):
        """ filter links by interval, end position, and possible bin conflict """
        # links_filtered = {}
        ctg2ctgs = {}
        n_too_small = 0
        _max_cov_gap = 0 #.3 * rlen
        _max_uncovered_end = max_insert_size - 2 * rlen
        n_kept = 0
        n_interval_len_filtered = 0
        n_max_cov_filtered = 0

        for (ctg1, ctg2), _pairs in _raw_links.items():
            if len(_pairs) < 4:
                n_too_small += 1
                continue
            ctg1_reads, ctg2_reads = zip(*_pairs)
            ctg1_len = ctg_lookup[ctg1].length
            ctg2_len = ctg_lookup[ctg2].length
            print(ctg1, ctg2)
            ctg1_covered_intervals, ctg1_nbases = get_cov_intervals(ctg1_reads, _max_cov_gap)
            ctg2_covered_intervals, ctg2_nbases = get_cov_intervals(ctg2_reads, _max_cov_gap)
            if ctg1_nbases and ctg2_nbases:
                ctg1_till_end, direction1 = _check_end_poz(ctg1_covered_intervals[0],
                                               ctg1_covered_intervals[1],
                                               ctg1_len)

                ctg2_till_end, direction2 = _check_end_poz(ctg2_covered_intervals[0],
                                               ctg2_covered_intervals[1],
                                               ctg2_len)
                if ctg1_till_end + ctg2_till_end < _max_uncovered_end:
                    bin1 = _get_bin_id(ctg1)
                    bin2 = _get_bin_id(ctg2)
                    interval1_len = _get_interval_len(ctg1_covered_intervals)
                    interval2_len = _get_interval_len(ctg2_covered_intervals)

                    # solve bin conflict by removing the weak links
                    if bin1 != bin2 or bin1 == -1 or bin2 == -1:
                        if max(interval1_len, interval2_len) < 500:
                            n_interval_len_filtered += 1
                            continue

                        if max(interval1_len/ctg1_len, interval2_len/ctg2_len) < 0.05:
                            n_max_cov_filtered += 1
                            continue

                    cov1 = _calc_cov(interval1_len, ctg1_nbases)
                    cov2 = _calc_cov(interval2_len, ctg2_nbases)
                    _update_ctg2ctgs(ctg2ctgs, ctg1, ctg2, interval1_len, cov1, bin2, direction1)
                    _update_ctg2ctgs(ctg2ctgs, ctg2, ctg1, interval2_len, cov2, bin1, direction2)
                    n_kept += 1

        print("after filter_links: ", n_kept)
        print("interval len filtered: ", n_interval_len_filtered)
        print("max cov perc filtered: ", n_max_cov_filtered)

        return ctg2ctgs

    def _pick_links_binned_one_end(_ctg:str, _ctg_bin, _linked_ctg:List[LinkedCtg], _todel: Set[Tuple[str, str]]):
        if _linked_ctg:
            picked = []
            unbinned = []
            diff_bin = []
            for il in _linked_ctg:
                if il.binid == _ctg_bin:
                    picked.append(il)
                elif il.binid > 0:
                    diff_bin.append(il)
                else:
                    unbinned.append(il)

            picked_bin = 0
            if not picked and unbinned:
                picked_bin = -1
                picked.append(max(unbinned, key=lambda x: x.interval_len * x.cov))

            if picked:
                if diff_bin:
                    for i in diff_bin:
                        _todel.update({(i.id, _ctg), (_ctg, i.id)})
                        print("=======deleted: (_ctg, i.id)")
                return picked, picked_bin
            else:
                if diff_bin:
                    picked.append(max(diff_bin, key=lambda x: x.interval_len * x.cov))
                    return picked, picked[0].binid

        else:
            return [], 0

    def _pick_links_unbinned_one_end(_ctg:str, _linked_ctg: List[LinkedCtg],
                                     _todel: Set[Tuple[str, str]]):
        if _linked_ctg:
            chosen = max(_linked_ctg, key=lambda x: x.interval_len * x.cov)
            if chosen.binid > 0:
                diff_bin = [il for il in _linked_ctg if 0 < il.binid != chosen.binid]
                if diff_bin:
                    for i in diff_bin:
                        _todel.update({(i.id, _ctg), (_ctg, i.id)})
                        print("=======deleted:", (i.id, _ctg))

            return [chosen], chosen.binid

        else:
            return [], 0

    def resolve_linked_cov(_ctg2ctgs: Dict[str, List[LinkedCtg]]):
        """ filter links by interval, end position, and possible bin conflict """
        _ctg_links = {}

        #=========================================== Testing
        # metabat_bin1 = _get_bin_id("NODE_35_length_245626_cov_16.939265")
        # metabat_bin5 = _get_bin_id("NODE_115_length_108006_cov_14.903786")

        for _ctg, linked_ctgs in _ctg2ctgs.items():
            _ctg_bin = _get_bin_id(_ctg)
            # ===========================================> Testing
            # if _ctg_bin == metabat_bin1:
            #     print("Found in bin1", _ctg)
            #     for il_ctg in linked_ctgs:
            #         print(il_ctg.id, il_ctg.cov, il_ctg.interval_len,il_ctg.binid)
            #
            # if _ctg_bin == metabat_bin5:
            #     print("Found in bin5", _ctg)
            #     for il_ctg in linked_ctgs:
            #         print(il_ctg.id, il_ctg.cov, il_ctg.interval_len,il_ctg.binid)
            # <=========================================== Testing

            if len(linked_ctgs) == 1:
                _ctg_links[(_ctg, linked_ctgs[0].id)] = linked_ctgs[0].interval_len
                # pair_key = tuple(sorted([_ctg, linked_ctgs[0].id]))
                # _ctg_links[pair_key] = linked_ctgs[0].nbases
                continue

            # ===========================================> Testing
            print("\n>>>",_ctg, _ctg_bin)
            for il_ctg in linked_ctgs:
                print(il_ctg.id, il_ctg.direction, il_ctg.cov, il_ctg.interval_len, il_ctg.binid)
            # <=========================================== Testing
            end_plus = []
            end_minus = []
            to_del = set()
            for il_ctg in linked_ctgs:
                if il_ctg.direction == 1:
                    end_plus.append(il_ctg)
                else:
                    end_minus.append(il_ctg)
            if _ctg_bin != -1:
                chosen_plus, bin_plus = _pick_links_binned_one_end(_ctg, _ctg_bin, end_plus, to_del)
                chosen_minus, bin_minus = _pick_links_binned_one_end(_ctg, _ctg_bin, end_minus, to_del)

                # ===========================================> Testing
                print("================= picking for binned: ")
                print([_i.id for _i in chosen_plus], bin_plus)
                print([_i.id for _i in chosen_minus], bin_minus,"\n")
                # ===========================================> Testing

                if 0 < bin_plus != bin_minus > 0:
                    # if the two ends are linked to 2 different bins and neither is _ctg_bin
                    chosen_plus = []
                    chosen_minus = []

            else:
                chosen_plus, bin_plus = _pick_links_unbinned_one_end(_ctg, end_plus, to_del)
                chosen_minus, bin_minus = _pick_links_unbinned_one_end(_ctg, end_minus, to_del)
                # ===========================================> Testing
                print("================= picking for unbinned: ")
                print([_i.id for _i in chosen_plus], bin_plus)
                print([_i.id for _i in chosen_minus], bin_minus, "\n")
                # <=========================================== Testing

            if chosen_plus:
                for il_plus in chosen_plus:
                    _ctg_links[tuple(sorted([_ctg, il_plus.id]))] = il_plus.interval_len
            if chosen_minus:
                for il_minus in chosen_minus:
                    _ctg_links[tuple(sorted([_ctg, il_minus.id]))] = il_minus.interval_len

            for idel in to_del:
                if idel in _ctg_links:
                    del _ctg_links[idel]

        # ===========================================> Testing
        print("============== Line 325, parsers, _ctg_links")
        for k, v in _ctg_links.items():
            print(k, v)
        # <=========================================== Testing

        return _ctg_links

    def _update_dict_val(k: Tuple[str, str], v: Tuple[ReadAlignment, ReadAlignment]) -> None:
        if k in raw_ctg_links:
            raw_ctg_links[k].append(v)
        else:
            raw_ctg_links[k] = [v]

    bam_handle = AlignmentFile(aln_file)
    raw_ctg_links = {}
    waiting_read_dict = {}
    for i_read in bam_handle.fetch():
        if not i_read.is_paired or i_read.is_supplementary:
            # or i_read.is_secondary \
            # if read.is_secondary or read.is_supplementary:
            continue

        i_rid = i_read.query_name
        if i_rid not in waiting_read_dict:
            if i_read.is_read1:
                waiting_read_dict[i_rid] = (i_read, None)
            elif i_read.is_read2:
                waiting_read_dict[i_rid] = (None, i_read)

        else:
            _imate = waiting_read_dict.pop(i_rid)
            i_mate = None
            if i_read.is_read1 and _imate[1]:
                i_mate = _imate[1]
            elif i_read.is_read2 and _imate[0]:
                i_mate = _imate[0]
            if not i_mate:
                continue

            if i_read.mapping_quality < 1 or i_mate.mapping_quality < 1:
                continue

            i_aln_len = i_read.query_alignment_length
            i_mate_len = i_mate.query_alignment_length

            if i_aln_len < aln_query_perc * i_read.query_length or \
                    i_mate_len < aln_query_perc * i_mate.query_length:
                continue

            i_ctg_id = i_read.reference_name
            i_mate_ctg_id = i_mate.reference_name

            if i_ctg_id in ctg_lookup and \
                    i_mate_ctg_id in ctg_lookup and \
                    i_ctg_id != i_mate_ctg_id:

                i_read_alignment = ReadAlignment(i_read.reference_start, i_read.reference_end)
                i_mate_alignment = ReadAlignment(i_mate.reference_start, i_mate.reference_end)
                _update_dict_val(tuple(sorted([i_ctg_id, i_mate_ctg_id])),
                                 (i_read_alignment, i_mate_alignment))
    bam_handle.close()
    ctg2ctg_links = filter_links(raw_ctg_links)
    filtered_links = resolve_linked_cov(ctg2ctg_links)

    return filtered_links
