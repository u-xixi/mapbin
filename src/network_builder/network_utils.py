from infomap import Infomap
from pysam import AlignmentFile, AlignedSegment
from typing import Dict, List, Tuple, TextIO
import time
import glob
from file_io import *
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
        print("[=== Binning layer ===] Processing binning result took {:.1f} secs".format(time_finish - time_begin))
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
                    # no need to sort the keys because the mate always shows up before the read
                    # thus always i_nbctg.num_id > i_mate_nbctg.num_id
                    _update_dict_val((i_nbctg.num_id, i_mate_nbctg.num_id), ij_weight, ctg_links)
        sam_handle.close()
        return ctg_links

    time_begin = time.time()
    print("[:::Message:::] Parsing alignment...")

    links = process_sam()
    n_edges = 0
    for (i, j), w in links.items():
        if w > 0.001:
            im_network.add_multilayer_intra_link(layer_idx, i, j, w)
            n_edges += 1
    time_finish = time.time()

    print("[=== Read pairing layer ===] Made {} edges based on read pairing.".format(n_edges))
    print("[=== Read pairing layer ===] Processing read pairing took {:.1f} secs".format(
        time_finish - time_begin))


def make_infomap_linked(args, im_network: Infomap, layer_idx:int, ctg_lookup) -> None:
    """
    does 3 steps:
    1. get read-to-barcode mapping
    2. get {barcode: (ctg_num_id, aligned_bases, ctg_length)}
    3. establish ctg-to-ctg edges.
    The edge weight is based on aligned_bases to total_bases
    The edge weight has to be > 0.001 (for now)
    """

    def _get_reads2bars(seq_file: str):
        seqfile = SeqFile(seq_file)
        seq_idf = seqfile.idf
        reads2bars = {}

        _header_spliter = lambda x: x[1:].split(" ")
        for iline in seqfile.file_content:
            if iline.startswith(seq_idf):
                iline_split = _header_spliter(iline)
                reads2bars[iline_split[0]] = iline_split[-1]
        return reads2bars

    def _update_ctg_dict(dict_:Dict[str, int], elem_key: str, elem_val: int)-> None:
        """when iter over each contig, a dict {bar: weight} is made.
        this call updates the dict"""
        if elem_key in dict_:
            dict_[elem_key] += elem_val
        else:
            dict_[elem_key] = elem_val

    def _update_bars2ctgs(b2c: Dict[str, List[Tuple[int, int, int]]], ctg_dict: Dict[str, int], ctg_id: str) -> None:
        """b2c: {bar: (ctg_num_id, mapped_bases, ctg_len)}
        ctg_dict: {bar: n_mapped_bases}"""
        ctg = ctg_lookup[ctg_id]
        ctg_num_id = ctg.num_id
        ctg_len = ctg.length
        for k, v in ctg_dict.items():
            if k in b2c:
                b2c[k].append((ctg_num_id, v, ctg_len))
            else:
                b2c[k] = [(ctg_num_id, v, ctg_len)]


    def get_bars2ctgs():
        b2c = {}  # bars2contigs, {barcode: [contigs]}
        reads2bars = {}
        if args.se_reads:
            reads2bars = _get_reads2bars(args.se_reads)
        if args.pe_reads:
            reads2bars.update(_get_reads2bars(args.pe_reads[0]))
        print("[:::Message:::] Processed {} reads with unique ids.".format(len(reads2bars)))

        sam_handle = AlignmentFile(args.aln_file, 'rb')
        last_ctg = None
        i_ctg_dict = {}

        for i_read in sam_handle.fetch():
            if i_read.is_unmapped or i_read.is_secondary or i_read.is_supplementary:
                continue

            i_ctg = i_read.reference_name
            if i_ctg not in ctg_lookup:
                continue

            if not last_ctg:
                last_ctg = i_ctg
                i_aln_len = i_read.query_alignment_length
                i_bar = reads2bars[i_read.query_name]
                i_ctg_dict[i_bar] = i_aln_len
            elif i_ctg != last_ctg:
                # conclude the last contig: update the b2c
                _update_bars2ctgs(b2c, i_ctg_dict, last_ctg)
                # set up the new contig
                i_ctg_dict = {}
                last_ctg = i_ctg
                i_aln_len = i_read.query_alignment_length
                i_bar = reads2bars[i_read.query_name]
                i_ctg_dict[i_bar] = i_aln_len
            else:
                i_aln_len = i_read.query_alignment_length
                i_bar = reads2bars[i_read.query_name]
                _update_ctg_dict(i_ctg_dict, i_bar, i_aln_len)

        sam_handle.close()

        return b2c

    def _make_ctg_link_single_bar(ctg_links: Dict[Tuple[int, int], float],
                                  ictgs: List[Tuple[int, int, int]]) -> None:
        nctgs = len(ictgs)
        for i in range(nctgs):
            for j in range(i, nctgs):
                ik = ictgs[i][0]
                jk = ictgs[j][0]
                k = None
                if ik < jk:
                    k = (ik, jk)
                else:
                    k = (jk, ik)

                v = (ictgs[i][1]+ ictgs[j][1])/(ictgs[i][2]+ ictgs[j][2])
                if k in ctg_links:
                    ctg_links[k] += v
                else:
                    ctg_links[k] = v

    def generate_ctg2ctg_edges(bars2ctgs_dict:Dict[str, List[Tuple[int, int, int]]],
                               min_bar_size: int) -> Dict[Tuple[int, int], float]:
        ctg_links = {}
        test_nbars = 0
        for ibar, ictgs in bars2ctgs_dict.items():
            if len(ictgs) >= min_bar_size:
                _make_ctg_link_single_bar(ctg_links, ictgs)
                test_nbars += 1

        print("[:::Dev Message:::] no. legit bars: {}".format(test_nbars))
        return ctg_links

    print("[=== Linked read layer ===] Parsing linked reads...")
    time_begin = time.time()

    bars2ctgs = get_bars2ctgs()
    print("[:::Message:::] processed reads with {} barcodes".format(len(bars2ctgs)))
    time_bars2ctgs = time.time()
    print("[:::Dev Message:::]bars2ctgs: {:.1f} secs".format(time_bars2ctgs - time_begin))
    # bars2ctgs elements:
    # BX:Z:A61B28C05D34: [(1, 0.00021912252873939136), (6, 0.0005723737009561289), (2428, 0.1162340178225494)]
    links = generate_ctg2ctg_edges(bars2ctgs, args.min_bar_size)
    nedges = 0
    for k, v in links.items():
        if v > 0.001:
            im_network.add_multilayer_intra_link(layer_idx, k[0], k[1], v)
            nedges += 1

    time_finish = time.time()
    print("[=== Linked read layer ===] made {} edges based on read barcodes".format(nedges))
    print("[=== Linked read layer ===] Processing read linkage took {:.1f} secs".format(time_finish - time_begin))


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


def generate_ctg_lookup(contig_fa_lookup: Dict[str, str], min_ctg_len: int):
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
            unbinned_short.append(iheader)
        else:
            formatted_id = format_ctg_idf(iheader)
            cur_nbcontig = NbContig("", ilen, i)  # changed node_id to int
            ctg_lookup[iheader] = cur_nbcontig  # changed the key to the fasta header
            im_network.set_name(i, formatted_id)
            # im_network.add_node(node_id=i,
            #                     name=formatted_id)

        i += 1
    print("[=== Parsing contigs ===] {} contigs with minimum length {} for binning".format(len(ctg_lookup), min_ctg_len))


    return ctg_lookup, im_network, unbinned_short




