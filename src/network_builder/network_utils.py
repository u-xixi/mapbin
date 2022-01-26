from infomap import Infomap
from Bio import SeqIO
from pysam import AlignmentFile, AlignedSegment
from typing import Dict, List, Tuple
import time
import sys
import glob
from file_io import *
from itertools import combinations
from collections import defaultdict
from parsers import *


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


format_ctg_idf = lambda x: "_".join(x.strip().split(" ")[0].split("_")[:2])

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


def make_infomap_binned(args, im_network:Infomap, layer_idx, ctg_lookup):

    def _get_bin_info_content(bin_fa):
        # returns a set
        print("[::Dev Message::] parsing \"{}\"".format(bin_fa))
        fa_content = _get_seq_file_content(bin_fa)
        fa_headers = [_parse_seq_name(i) for i in fa_content if i.startswith(">")]  # headers don't have ">"
        fa_headers = [i for i in fa_headers if i in ctg_lookup]
        n_edges = 0

        if len(fa_headers) > 1:

            for a, b in combinations(fa_headers, 2):
                ctg_a = ctg_lookup[a]
                ctg_b = ctg_lookup[b]
                # print(ctg_a)
                a_key = ctg_a.num_id
                b_key = ctg_b.num_id
                a2b_val = args.min_ctg_len / ctg_a.length
                # b2a_key = " ".join([ctg_b.num_id, ctg_a.num_id])
                b2a_val = args.min_ctg_len / ctg_b.length
                links = [((layer_idx, a_key), (layer_idx, b_key), a2b_val),
                         ((layer_idx, b_key), (layer_idx, a_key), b2a_val)]
                # links = [((layer_idx, ctg_a), (layer_idx, ctg_b)),
                #          ((layer_idx, ctg_b), (layer_idx, ctg_a))]
                n_edges += 2
                im_network.add_multilayer_links(links)
            print("[::Dev Message::] added {} edges".format(n_edges) )
            return 1
        else:
            return 0
        # return im_edges

    def _get_seq_file_content(file_path):

        seq_file = SeqFile(file_path)
        if seq_file.seqformat == "fasta":
            return seq_file.file_content
        else:
            sys.exit("Invalid Fasta file: {}".format(file_path))

    # def _get_fasta_files_in_dir(dir, suffix):
    #     # dir must ends up with "/"
    #     return glob.glob(dir + "*." + suffix)

    _parse_seq_name = lambda x: "_".join(x[1:-1].split("_")[:2])

    for i_bin_dir in args.bin_dir:
        bin_fa_files = glob.glob(i_bin_dir + "*." + args.bin_suffix)
        print("[::Message::] Parsing binning result directory \"{}\".".format(i_bin_dir))
        if not len(bin_fa_files):
            print("[::Warning::] No contigs detected in \"{}\". Directory skipped".format(i_bin_dir))
            continue
        # i_edges = {}

        time_begin = time.time()
        print("[::Message::] Parsing binning results...")

        sig = 0
        for bin_fa in bin_fa_files:
            sig += _get_bin_info_content(bin_fa)

        if sig == 0:
            print("[::Warning::] No useful binning result in \"{}\". Directory skipped".format(i_bin_dir))

        # update_dict_val(im_edges_content, i_edges)
        time_finish = time.time()
        print("[=== Binning layer ===] parsing binning result takes {:.1f} secs".format(time_finish - time_begin))


def make_infomap_pairing(args, im_network: Infomap, layer_idx, ctg_lookup):

    def _update_edge_weight(k, v):
        nonlocal ctg_net
        if k in ctg_net:
            ctg_net[k] += v
        else:
            ctg_net[k] = v

    def _get_ctg_info(ctg_raw_str: str) -> (str, int):
        # formatted_id = "_".join(ctg_raw_str.split("_")[:2])
        formatted_id = format_ctg_idf(ctg_raw_str)
        if formatted_id in ctg_lookup:
            nbctg = ctg_lookup[formatted_id]
            return nbctg.num_id, nbctg.length
        return None, None

    def _conclude_last_pair(r1_: AlignedSegment, r2_: AlignedSegment) -> None:
        # get the list of ctg that r1 and r2 share
        # this function changes ctg_net
        # nonlocal ctg_net
        nonlocal test_counter
        r1_ctg_ = r1_.reference_name
        r2_ctg_ = r2_.reference_name

        if r1_ctg_ != r2_ctg_:
            test_counter += 1
            # print("[:::Dev Message:::] find pair with diff parents +1")

            r1_ctg_num_id, r1_ctg_len = _get_ctg_info(r1_ctg_)
            r2_ctg_num_id, r2_ctg_len = _get_ctg_info(r2_ctg_)

            if r1_ctg_num_id and r2_ctg_num_id:
                ctg_pair_key12 = (r1_ctg_num_id, r2_ctg_num_id)
                ctg_pair_key21 = (r2_ctg_num_id, r1_ctg_num_id)
                edge_wt12 = r1.reference_length/r1_ctg_len
                edge_wt21 = r2.reference_length/r2_ctg_len
                _update_edge_weight(ctg_pair_key12, edge_wt12)
                _update_edge_weight(ctg_pair_key21, edge_wt21)

            # this means the pair is gonna be kept
            # in this case, first, format the r1_ctg and r2_ctg strings



    # def write_pairing_layer(ctg_net_: Dict):
    #     ctg_net_as_list = [" ".join([" ".join(k), "{:.5f}".format(v)]) for k, v in ctg_net_.items()]
    #     ctg_net_as_list.sort()
    #     del ctg_net_
    #     return "\n".join(ctg_net_as_list) + "\n"

    time_begin = time.time()
    print("[::Message::] Parsing alignment...")

    alnfile = AlignmentFile(args.aln_file, 'r')
    ctg_net = {}

    waiting_read_dict = defaultdict(lambda: [None, None])
    # r1 = None
    # r2 = None
    test_counter = 0

    for read in alnfile.fetch():

        if not read.is_proper_pair or read.is_secondary or read.is_supplementary:
            # if read.is_secondary or read.is_supplementary:
            continue

        read_id = read.query_name
        if read_id not in waiting_read_dict:
            if read.is_read1:
                waiting_read_dict[read_id][0] = read
            else:
                waiting_read_dict[read_id][1] = read

        else:
            if read.is_read1:
                r1 = read
                r2 = waiting_read_dict[read_id][1]
                _conclude_last_pair(r1, r2)
                if waiting_read_dict[read_id][0]:
                    print("[:::Dev Warning:::] pair abnormal: 2x R1 of read \"{}\"".format(read_id))
            else:
                r2 = read
                r1 = waiting_read_dict[read_id][0]
                _conclude_last_pair(r1, r2)
                if waiting_read_dict[read_id][1]:
                    print("[:::Dev Warning:::] pair abnormal: 2x R2 of read \"{}\"".format(read_id))
            del waiting_read_dict[read_id]
            # test_counter += 1


        # if read.is_read1:
        #     if r1 and r2:
        #         _conclude_last_pair(r1, r2, test_counter)
        #         # test_counter += 1
        #
        #     if read.is_unmapped or read.mate_is_unmapped:
        #         r1 = None
        #     else:
        #         r1 = read
        # elif read.is_read2:
        #     if not r1:
        #         r2 = None
        #     else:
        #         r2 = read

    print("[::Dev Message::]Line 185: test_counter = {} ".format(test_counter))
    print("[::Dev Message::]waiting read dict size {}".format(len(waiting_read_dict)))
    alnfile.close()

    time_checkpoint1 = time.time()
    print("[=== Read pairing layer ===] parsing alignment takes {:.1f} secs".format(time_checkpoint1 - time_begin))
    print("[::Dev Message::]Line 189: size of ctg_net = {} ".format(len(ctg_net)))

    if ctg_net:
        for k, v in ctg_net.items():
            multilayer_link_source = (layer_idx, k[0])
            multilayer_link_tgt = (layer_idx, k[1])
            im_network.add_multilayer_link(multilayer_link_source, multilayer_link_tgt, v)
        n_edges = len(ctg_net)
        print("[::Dev Message::] added {} edges".format(n_edges))


    time_checkpoint2 = time.time()
    print("[=== Read pairing layer ===] Build infomap layer takes {:.1f} secs".format(time_checkpoint2- time_checkpoint1))



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
                im_network.add_multilayer_link(source_multilayer_node=(layer_idx, ikey[0]),
                                               target_multilayer_node=(layer_idx, ikey[1]),
                                               weight=ival)
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


def make_infomap_assembly(args, im_network: Infomap, layer_idx: int, all_contigs: Set[str], ctg_lookup: Dict[str, NbContig]):
    # get paths: {contig_num_id: [(path_left_end, path_right_end)]}
    #  Dict[int, List[Tuple]]
    time_begin = time.time()
    if args.assembler == "spades":
        path_segs = parse_spades_contig_paths(all_contigs, ctg_lookup, args.ctg_paths)
        contigs = list(path_segs.keys())
        contigs.sort()
        print(len(contigs))
        counter = 0
        for i in range(len(contigs)):
            for j in range(i+1, len(contigs)):
                ic = contigs[i]
                jc = contigs[j]
                isegs = path_segs[ic]
                jsegs = path_segs[jc]
                if isegs.intersection(jsegs):
                    im_network.add_multilayer_link(source_multilayer_node=(layer_idx, ic),
                                                   target_multilayer_node=(layer_idx, jc))
                    counter +=1
        print("[=== Assembly layer ===] linked {} contigs".format(counter))
    time_end = time.time()

    print("[:::Dev Message:::]assembly processing: {:.1f} secs".format(time_end - time_begin))







# tqdm python for printing out the progress
def generate_ctg_lookup(ctg_fa, min_ctg_len):
    # function: generate ctg lookup dict
    #           and initialize im by initializing the network and adding nodes
    # change: node_id is no longer a str, but an int
    ctg_lookup = {}
    unbinned = []
    i = 1
    im_network = Infomap("-i multilayer --two-level")

    for seq_record in SeqIO.parse(ctg_fa, "fasta"):
        seq_len = len(seq_record)
        if seq_len < min_ctg_len:
           unbinned.append(seq_record)
        else:
            formatted_id = "_".join(seq_record.id.split("_")[:2])
            cur_nbcontig = NbContig(formatted_id, seq_len, i)  # changed node_id to int
            ctg_lookup[cur_nbcontig.id] = cur_nbcontig
            im_network.add_node(node_id=i,
                                name=formatted_id)

        i += 1
    return ctg_lookup, unbinned, im_network


def generate_ctg_lookup_alt(contig_fa_lookup: Dict[str, str], min_ctg_len: int):
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
    :return: im_network: the infomap network built by this func
    """
    ctg_lookup = {}
    unbinned = []
    i = 1
    im_network = Infomap("-s 1")

    for iheader, iseq in contig_fa_lookup.items():
        ilen = len(iseq)

        if ilen < min_ctg_len:
            unbinned.append("\n".join([">" + iheader, iseq]))
        else:
            formatted_id = "_".join(iheader.split("_")[:2])
            cur_nbcontig = NbContig(formatted_id, ilen, i)  # changed node_id to int
            ctg_lookup[iheader] = cur_nbcontig  # changed the key to the fasta header
            im_network.add_node(node_id=i,
                                name=formatted_id)

        i += 1
    return ctg_lookup, unbinned, im_network


def run_infomap(im_network:Infomap):
    im_network.run()


