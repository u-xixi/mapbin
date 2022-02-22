from infomap import Infomap
from Bio import SeqIO
from pysam import AlignmentFile, AlignedSegment
from typing import Dict, List, Tuple, TextIO
import time
import sys
import glob
from file_io import *
from itertools import combinations
from collections import defaultdict
from parsers import *
from random import shuffle, seed


class NbContig():
    # contig is the basic unit of this network
    # use this class as the basis
    # it might help reduce the overhead
    # and allows modifying or adding info to the contig
    __slots__ = ['id', 'length', 'num_id']
    def __init__(self, id:str, length=NotImplemented, num_id=NotImplemented):
        """
        :param id: NODE_1
        :param length: contig length
        :param num_id: numeric id for informap
        """
        self.id = id
        self.length = length
        self.num_id = num_id


# format_ctg_idf = lambda x: "_".join(x.strip().split(" ")[0].split("_")[:2])
format_ctg_idf = lambda x: x.strip().split(" ")[0]

# def make_infomap_nodes(ctg_lookup: Dict[str, NbContig]):
#     "multilayer network builder using Infomap API"
#     im_network = Infomap("-i ")
#     for i_nbc in ctg_lookup.items():
#         im_network.add_node(i_nbc)
    # with open(contig_fa, 'r') as fr:
    #     for iline in fr
#     add nodes from the contig list
#     input: contig.fasta
#     output: im_network


def make_complete_graph(nodes: List[NbContig], im_network, layer_idx) -> str:
    if len(nodes) == 1:
        return ""

    else:
        obuffer = []
        link_count = 0
        for i in range(len(nodes)):
            for j in range(i+1, len(nodes)):
                inode = nodes[i]
                jnode = nodes[j]
                im_network.add_multilayer_intra_link(layer_idx, inode.num_id,
                                                     jnode.num_id, weight=1)
                link_count += 1
                obuffer.append("{0} {1} {2} {3}".format(layer_idx, inode.num_id, jnode.num_id, 1))
        print("[:::Message:::] made complete graph with {} edges".format(link_count))

        return "\n".join(obuffer)


def make_infomap_binned(args, im_network: Infomap, layer_idx,
                        ctg_lookup: Dict[str, NbContig]) -> str:

    def _subbin_handler(subbin: List[NbContig], obuffer:List[str]) -> NbContig:
        """
        decide which contig is the center (median) that connects to all the rest of contigs
        which contig is the representative (biggest) that serves as the node connected to the other subbins.
        this handler makes star
        :return:
        """
        if len(subbin) == 1:
            return subbin[0]

        ctg_lens = [i.length for i in subbin]
        ctg_lens.sort()
        ctg_median = ctg_lens[len(ctg_lens) // 2 - 1]
        ctg_max = ctg_lens[-1]

        center_ctg = None
        rep_ctg = None
        circle_ctg = []

        for i in subbin:
            ilen = i.length
            if ilen == ctg_median:
                center_ctg = i
            elif ilen == ctg_max:
                rep_ctg = i
                circle_ctg.append(rep_ctg)
            else:
                circle_ctg.append(i)

        for i_ctg in circle_ctg:
            im_network.add_multilayer_intra_link(layer_idx, i_ctg.num_id, center_ctg.num_id, weight=1)
            # to be solved: weight between nodes
            obuffer.append("{0} {1} {2} {3}".format(layer_idx, i_ctg.num_id, center_ctg.num_id, 1))
        return rep_ctg

    def make_network_bin_fa(bin_fa: str, s: int) -> str:
        """

        :param bin_fa: fasta file represent one bin
        :param s: number of sub-bins
        :return:
        """
        # step 1: decide the bin-dividing plan
        # step 2: fill in the edges


        # step 1
        bin_headers = set(parse_contig_fa(bin_fa).keys())
        bin_ctg = [ctg_lookup[i] for i in bin_headers if i in ctg_lookup]
        # return make_complete_graph(nodes=bin_ctg, im_network=im_network, layer_idx=layer_idx)

        if len(bin_ctg) <= s:
            print("[:::Message:::] make complete graph for {}".format(bin_fa))
            return make_complete_graph(nodes=bin_ctg, im_network=im_network, layer_idx=layer_idx)

        else:
            print("[:::Message:::] make density-reduced graph for {}".format(bin_fa))
            bin_obuffer = []
            seed(12345)
            shuffle(bin_ctg)

            N = len(bin_ctg)
            Lmin = N//s
            Lmax = Lmin + 1
            all_L = [Lmax] * (N%s)
            all_L.extend([Lmin]*(s - N%s))
            # print("[:::Dev Message:::]network_utils, L133, no. nodes: {}; sizes: {}".format(N, all_L))
            # print("[:::Dev Message:::]network_utils, L133, no. nodes: {}\n {}".format(bin_fa, N))

            # step 2
            l0 = 0
            reps = []
            for l1 in all_L:
                i_subbin = bin_ctg[l0:l0+l1]
                irep = _subbin_handler(i_subbin, obuffer)
                reps.append(irep)
                l0 += l1

            buffer_text_p2 = make_complete_graph(reps, im_network=im_network,
                                layer_idx=layer_idx)
            bin_obuffer.append(buffer_text_p2)

            return "\n".join(bin_obuffer)


# def make_infomap_binned(args, im_network:Infomap, layer_idx, ctg_lookup, out_fhandle: TextIO)-> None:
    obuffer = []
    for i_bin_dir in args.bin_dir:
        bin_fa_files = glob.glob(i_bin_dir + "*." + args.bin_suffix)
        print("[=== Binning layer ===] Parsing binning result directory \"{}\".".format(i_bin_dir))
        if not len(bin_fa_files):
            print("[::Warning::] No contigs detected in \"{}\". Directory skipped".format(i_bin_dir))
            continue

        time_begin = time.time()

        for i_fa in bin_fa_files:
            i_buffer_text = make_network_bin_fa(bin_fa=i_fa, s=args.s)
            obuffer.append(i_buffer_text)

        # update_dict_val(im_edges_content, i_edges)
        time_finish = time.time()
        print("[=== Binning layer ===] Processing binning result takes {:.1f} secs".format(time_finish - time_begin))
        return "\n".join(obuffer)


def make_infomap_pairing(args, im_network: Infomap, layer_idx:int, ctg_lookup:Dict[str, NbContig]):
    def _update_dict_val(k: Tuple[int, int], v: int, ctg_links: Dict[Tuple[int, int], int]) -> None:
        if k in ctg_links:
            ctg_links[k] += v
        else:
            ctg_links[k] = v

    def process_sam() -> Dict[Tuple[int, int], int]:
        sam_handle = AlignmentFile(args.aln_file, 'rb')
        ctg_links = {}
        waiting_read_dict = {}
        for i_read in sam_handle.fetch():
            if not i_read.is_paired or i_read.is_secondary or i_read.is_supplementary:
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

                i_ctg_id = i_read.reference_name
                i_mate_ctg_id = i_mate.reference_name
                if i_ctg_id in ctg_lookup and \
                        i_mate_ctg_id in ctg_lookup and \
                        i_ctg_id != i_mate_ctg_id:
                    i_nbctg = ctg_lookup[i_ctg_id]
                    i_mate_nbctg = ctg_lookup[i_mate_ctg_id]
                    ij_weight = (i_read.query_alignment_length +
                                 i_mate.query_alignment_length)/(i_nbctg.length + i_mate_nbctg.length)
                    _update_dict_val((i_nbctg.num_id, i_mate_nbctg.num_id), ij_weight, ctg_links)
        sam_handle.close()
        return ctg_links

    time_begin = time.time()
    print("[:::Message:::] Parsing alignment...")

    links = process_sam()
    n_edges = 0
    for (i, j), w in links.items():
        if w > 0.00001:
            im_network.add_multilayer_intra_link(layer_idx, i, j, w)
            n_edges += 1
    time_finish = time.time()

    print("[=== Read pairing layer ===] Made {} edges based on read pairing.".format(n_edges))
    print("[=== Read pairing layer ===] Processing read pairing takes {:.1f} secs".format(
        time_finish - time_begin))


def make_infomap_linked(args, im_network: Infomap, layer_idx:int, ctg_lookup):

    def _get_reads2bars():
        seqfile = SeqFile(seq_file_path)
        seq_idf = seqfile.idf
        reads2bars = {}

        _header_spliter = lambda x: x[1:].split(" ")
        for iline in seqfile.file_content:
            if iline.startswith(seq_idf):
                iline_split = _header_spliter(iline)
                reads2bars[iline_split[0]] = iline_split[-1]
        return reads2bars

    _update_ctg_dict = lambda dict_, elem_key, elem_val: dict_[elem_key] + elem_val if elem_key in dict_ else elem_val

    def _update_bars2ctgs(base_dict, update_dict):
        # update_elem = lambda elem_k, elem_v: base_dict.setdefault(elem_k, []).append(elem_v)
        #  the barcode has a list of (contigs, mapped_base_perc)
        for key, val in update_dict.items():
            if key in base_dict:
                base_dict[key].append(val)
            else:
                base_dict[key] = [val]
            # update_elem(key, val)

    def _bars2ctgs_fltr(bars2ctgs_dict):
        bars2rm = []
        for bar, ctgs in bars2ctgs_dict.items():
            if len(ctgs) < 2:
                bars2rm.append(bar)

        for bar in bars2rm:
            del bars2ctgs_dict[bar]
        # return bars2ctgs

    def _get_bars2ctgs(alignment_file:str):
        # nonlocal aln_file
        # nonlocal ctg_lookup
        bars2ctgs = {}  # {barcode: [contigs]}
        # test_counter = 0
        reads2bars = _get_reads2bars()
        # print("[:::Dev Message:::] Line 278, reads2bars size: {}".format(len(reads2bars)))
        # checkpoint: successful till here

        aln_handle = AlignmentFile(alignment_file, 'rb')
        # print("[:::Dev Message:::] Line 281, no. alignment ref: {}".format(len(aln_handle.references)))
        # print("[:::Dev Message:::] Line 282, ctg lookup has {} items".format(len(ctg_lookup)))

        # check if the formatting of contig names works: debugged
        # test_ctg = sample(aln_handle.references, 1)[0]
        # test_ctg_formatted = format_ctg_idf(test_ctg)
        # test_ctg_in_lookup = sample(ctg_lookup.keys(), 1)[0]
        # print("[:::Dev Message:::] Line 291, ctg from aln file: {}\nctg from ctg_lookup: {}"
        #       .format(test_ctg_formatted, test_ctg_in_lookup))

        for i_ctg in aln_handle.references:
            i_formatted_ctg = format_ctg_idf(i_ctg)

            if i_formatted_ctg in ctg_lookup:
                # tho it is an inner function, it can access ctg_lookup
                i_ctg_len = ctg_lookup[i_formatted_ctg].length
                i_ctg_numid = ctg_lookup[i_formatted_ctg].num_id
                ctg_dict = {}  # {barcode: aligned_bases}

                for iread in aln_handle.fetch(i_ctg):
                    if iread.is_unmapped or iread.is_secondary or iread.is_supplementary:
                        continue

                    cur_bar = reads2bars[iread.query_name]
                    if cur_bar in ctg_dict:
                        # print("[:::Dev Message:::] Line 307, in ctg_dict:\n{}: {}".format(cur_bar, ctg_dict[cur_bar]))
                        ctg_dict[cur_bar] += iread.reference_length
                    else:
                        ctg_dict[cur_bar] = iread.reference_length

                if ctg_dict:
                    # test_counter += 1  ## debugged
                    # now the dict is being "harvested"
                    # ctg_dict: {bar: (contig_name, alignment_weighted_by_contig_size)}
                    ctg_dict = {k: (i_ctg_numid, v / i_ctg_len) for k, v in ctg_dict.items()}
                    _update_bars2ctgs(bars2ctgs, ctg_dict)

        aln_handle.close()
        # print("[:::Dev Message:::]Line 312, test_counter = {}".format(test_counter))

        _bars2ctgs_fltr(bars2ctgs)
        return bars2ctgs

    def _update_ctg_network(ctg_network_dict: Dict, key_ ,val_):
    # def _update_ctg_network(key_, val_):
        if key_ in ctg_network_dict:
            ctg_network_dict[key_] += val_
        else:
            ctg_network_dict[key_] = val_

    def make_barcode_network(bars2ctgs_dict:Dict[str, List[Tuple[int, float]]]):
        ctg_network = {}  # (contig_x, contig_y) as key

        for bar, contig_list in bars2ctgs.items():
            # contigs = [_parse_contig_name(i.id) for i in contig_list]
            for combo in combinations(contig_list, 2):  # permutations returns both (a, b) and (b,a)
                key01 = (combo[0][0], combo[1][0])
                val01 = combo[0][1]
                key10 = (combo[1][0], combo[0][0])
                val10 = combo[1][1]

                _update_ctg_network(ctg_network, key01, val01)
                _update_ctg_network(ctg_network, key10, val10)

        n_edges = 0
        for ikey, ival in ctg_network.items():
            if ival > 0.001:
                n_edges += 1
                im_network.add_multilayer_intra_link(layer_idx, ikey[0], ikey[1], weight=ival)
        print("[:::Dev Message:::] added {} edges".format(n_edges))

    aln_file = args.aln_file
    seq_file_path = ""
    if args.se_reads:
        # need fixing: can only accept se reads atm
        seq_file_path = args.se_reads
    elif args.pe_reads:
        seq_file_path = args.pe_reads[0]

    # layer_id = str(layer_idx)
    # im_layer_content = set()
    time_begin = time.time()
    bars2ctgs = _get_bars2ctgs(aln_file)

    # for sanity check
    print("[:::Dev Message:::]bars2ctgs size: {}".format(len(bars2ctgs)))
    # sanity check passed
    # bars2ctgs elements:
    # BX:Z:A61B28C05D34: [(1, 0.00021912252873939136), (6, 0.0005723737009561289), (2428, 0.1162340178225494)]
    # icounter = 0
    # for k, v in bars2ctgs.items():
    #     print("{}: {}".format(k, v))
    #     icounter += 1
    #     if icounter > 2:
    #         break

    time_bars2ctgs = time.time()
    print("[:::Dev Message:::]bars2ctgs: {:.1f} secs".format(time_bars2ctgs - time_begin))

    make_barcode_network(bars2ctgs)
    time_ctg_net = time.time()
    print("[:::Dev Message:::]ctg network building: {:.1f} secs".format(time_ctg_net - time_bars2ctgs))


def make_infomap_assembly(args, im_network: Infomap,
                          layer_idx: int, all_contigs: Set[str],
                          ctg_lookup: Dict[str, NbContig]):
    # get paths: {contig_num_id: [(path_left_end, path_right_end)]}
    #  Dict[int, List[Tuple]]
    time_begin = time.time()
    # obuffer = []

    if args.assembler == "spades":
        path_segs = parse_spades_contig_paths(all_contigs, ctg_lookup, args.ctg_paths)
        contigs = list(path_segs.keys())
        contigs.sort()

        counter = 0
        for i in range(len(contigs)):
            for j in range(i+1, len(contigs)):
                ic = contigs[i]
                jc = contigs[j]
                isegs = path_segs[ic]
                jsegs = path_segs[jc]
                if isegs.intersection(jsegs):
                    im_network.add_multilayer_intra_link(layer_idx, ic, jc)
                    # obuffer.append("{0} {1} {2}".format(layer_idx, ic, jc))
                    counter +=1
        print("[=== Assembly layer ===] made {} links with {} contigs".format(counter, len(contigs)))
    time_end = time.time()

    print("[=== Assembly layer ===] Processing assembly took {:.1f} secs".format(time_end - time_begin))
    # return "\n".join(obuffer)


def generate_ctg_lookup(contig_fa_lookup: Dict[str, str], min_ctg_len: int, unbinned_short_file:str):
    # alternative func for generate ctg lookup
    # function: generate ctg lookup dict
    #           and initialize im by initializing the network and adding nodes
    # change: node_id is no longer a str, but an int
    """
    :param contig_fa_lookup: {fasta_header: seq_len}
    :param min_ctg_len:
    :return: contig_fa_lookup.keys(): all contig names as a set
    :return: ctg_lookup: Dict[str, NbContig]
            contigs that passed the min contig length check, formatted for use in infomap
    :return: unbinned: [fasta record as string]
    :return: im_network is no longer generated here.
    """
    im_network = Infomap("-s 1")
    ctg_lookup = {}
    unbinned_short = []
    i = 1

    for iheader, iseq in contig_fa_lookup.items():
        ilen = len(iseq)

        if ilen < min_ctg_len:
            unbinned_short.append("\n".join([">" + iheader, iseq]))
        else:
            formatted_id = format_ctg_idf(iheader)
            cur_nbcontig = NbContig("", ilen, i)  # changed node_id to int
            ctg_lookup[iheader] = cur_nbcontig  # changed the key to the fasta header
            im_network.set_name(i, formatted_id)
            # im_network.add_node(node_id=i,
            #                     name=formatted_id)

        i += 1
    print("[=== Parsing contigs ===] {} contigs with minimum length {} for binning".format(len(ctg_lookup), min_ctg_len))
    print("[=== Parsing contigs ===] writing out contigs below minimum length to file: {}" .format(unbinned_short_file))
    with open(unbinned_short_file, 'w') as fw:
        fw.write("\n".join(unbinned_short) + "\n")

    return ctg_lookup, im_network




