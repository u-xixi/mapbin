from itertools import combinations
from collections import deque
from file_io import *
# from Bio import SeqIO
import time
import pysam


def _get_reads2bars(seq_file_path):
    seqfile = SeqFile(seq_file_path)
    file_content = deque(seqfile.file_content)
    n_non_header = seqfile.step - 1
    _header_spliter = lambda x: x[1:].split(" ")
    reads2bars = {}
    while file_content:
        curline = file_content.popleft()
        _line_split = _header_spliter(curline)

        reads2bars[_line_split[0]] = _line_split[-1]
        for i in range(n_non_header):
            file_content.popleft()
    return reads2bars


# def get_bars2contigs(alnfile, seq_file_path, min_contig_len):
#     refs_and_lens = [(ctg, ctg_len) for (ctg, ctg_len) \
#                            in zip(alnfile.references, alnfile.lengths)\
#                            if ctg_len > min_contig_len]
#     reads2bars = _get_reads2bars(seq_file_path)
#
#     bars2contigs = {}  # {barcode: [list_of_NbNetContigs_aligned_by_this_barcode]}
#     _update_ref_dict = lambda dict_, key, val: dict_[key] + val if key in dict_ else val
#     for ctg, ctg_len in refs_and_lens:
#         ctg_dict = {}  # {bar: no_aligned_bases}, scope: one contig
#         for r in alnfile.fetch(ctg):
#             # r is a read aligned to ref
#             cur_bar = reads2bars[ctg.query_name]
#             # a contig might have reads that share the barcode, so the aligned bases accumulate
#             ctg_dict[cur_bar] = _update_ref_dict(ctg_dict, cur_bar, r.reference_length)
#         if ctg_dict:
#             # now the dict is being "harvested"
#             # ctg_dict: {bar: (contig_name, alignment_weighted_by_contig_size)}
#             ctg_dict = {k: (ctg, ctg_len/v) for k,v in ctg_dict.items()}
#             _update_list_dict(bars2contigs, ctg_dict)
#     return bars2contigs


def get_bars2ctgs(aln_file_path, seq_file_path, min_contig_len, CONTIG_LOOKUP, min_bar_size=2):
    alnfile = pysam.AlignmentFile(aln_file_path, 'rb')
    refs_and_lens = [(ctg, ctg_len) for (ctg, ctg_len) \
                           in zip(alnfile.references, alnfile.lengths)\
                           if ctg_len > min_contig_len]
    reads2bars = _get_reads2bars(seq_file_path)

    bars2ctgs = {}  # {barcode: [list_of_NbNetctgs_aligned_by_this_barcode]}
    _update_ctg_dict = lambda dict_, key, val: dict_[key] + val if key in dict_ else val
    _get_ctg_id = lambda x: CONTIG_LOOKUP["_".join(x.split("_")[:2])].num_id

    for ctg, ctg_len in refs_and_lens:
        ctg_dict = {}  # {bar: no_aligned_bases}, scope: one contig
        ctg_id = _get_ctg_id(ctg)

        for r in alnfile.fetch(ctg):
            # r is a read aligned to ref
            cur_bar = reads2bars[r.query_name]

            # a contig might have reads that share the barcode, so the aligned bases accumulate
            ctg_dict[cur_bar] = _update_ctg_dict(ctg_dict, cur_bar, r.reference_length)
        if ctg_dict:
            # now the dict is being "harvested"
            # ctg_dict: {bar: (contig_name, alignment_weighted_by_contig_size)}
            ctg_dict = {k: (ctg_id, v/ctg_len) for k,v in ctg_dict.items()}
        _update_bars2ctgs(bars2ctgs, ctg_dict)
    alnfile.close()
    _bars2ctgs_fltr(bars2ctgs, min_bar_size)
    return bars2ctgs


def _update_bars2ctgs(base_dict, update_dict):
    update_elem = lambda k, v: base_dict.setdefault(key, []).append(val)
    #  the barcode has a list of (contigs, mapped_base_perc)
    for key, val in update_dict.items():
        update_elem(key, val)


def _bars2ctgs_fltr(bars2ctgs, min_bar_size):
    bars2rm = []
    for bars, ctgs in bars2ctgs.items():
        if len(ctgs) < min_bar_size:
            bars2rm.append(bars)
    for bar in bars2rm:
        del bars2ctgs[bar]
    return bars2ctgs


def make_barcode_network(bars2ctgs):
    ctg_network = {}  # (contig_x, contig_y) as key
    def update_ctg_network(k, v):
        if k in ctg_network:
            ctg_network[k] += v
        else:
            ctg_network[k] = v
    # dict_val_updater = lambda key, val: ctg_network[key] + val if key in ctg_network else val

    for bar, contig_list in bars2ctgs.items():
        # contigs = [_parse_contig_name(i.id) for i in contig_list]
        for combo in combinations(contig_list, 2):  # permutations returns both (a, b) and (b,a)
            key1 = (combo[0][0], combo[1][0])
            val1 = combo[0][1]
            key2 = (combo[1][0], combo[0][0])
            val2 = combo[1][1]
            update_ctg_network(key1, val1)
            update_ctg_network(key2, val2)
    return ctg_network


def make_barcode_net(args, CONTIG_LOOKUP, layer_id):
    # required arg attributes:
    aln_file = args.aln_file
    seq_file_path = ""
    if args.se_reads:
        # need fixing: can only accept se reads atm
        seq_file_path = args.se_reads
    min_ctg_len = args.min_ctg_len
    min_bar_size = args.min_bar_size

    layer_id = str(layer_id)
    im_layer_content = set()
    time_begin = time.time()
    bars2ctgs = get_bars2ctgs(aln_file, seq_file_path, min_ctg_len, CONTIG_LOOKUP, min_bar_size)
    time_bars2ctgs = time.time()
    print("[bars2ctgs] {:.1f} secs".format(time_bars2ctgs - time_begin))

    ctg_network = make_barcode_network(bars2ctgs)
    time_ctg_net = time.time()
    print("[ctg network building] {:.1f} secs".format(time_ctg_net - time_bars2ctgs))

    for k in ctg_network:  # k is a contig pair ('0', '27')
        f1 = " ".join(k)  # layer_id, node1_id, node2_id, weight
        f2 = "{:5f}".format(ctg_network[k])
        line = " ".join([layer_id, f1, f2])
        im_layer_content.add(line)
    time_layer_set = time.time()
    print("[make layer set] {:.1f} secs".format(time_layer_set - time_ctg_net))

    im_layer_content = list(im_layer_content)
    im_layer_content.sort()

    return "\n".join(im_layer_content) + "\n"


#
# class NbNetContig():
#     # contig is the basic unit of this network
#     # use this class as the basis
#     # it might help reduce the overhead
#     # and allows modifying or adding info to the contig
#     __slots__ = ('id', 'weight')
#     def __init__(self, id, weight=NotImplemented):
#         self.id = id
#         self.weight = weight
