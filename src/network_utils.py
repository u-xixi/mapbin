from infomap import Infomap
from typing import Dict, List, Tuple, Set
import time
from parsers import *
from assembly_graphs import *
from math import sqrt,ceil, log2



# format_ctg_idf = lambda x: "_".join(x.strip().split(" ")[0].split("_")[:2])
format_ctg_idf = lambda x: x.strip().split(" ")[0]


def make_complete_graph(nodes: List[NbContig], im_network, layer_idx) -> None:
    if len(nodes) > 1:
        link_count = 0
        for i in range(len(nodes)):
            for j in range(i+1, len(nodes)):
                inode = nodes[i]
                jnode = nodes[j]
                im_network.add_multilayer_intra_link(layer_idx, inode.num_id,
                                                     jnode.num_id, weight=1)
                link_count += 1
        print("[=== Info ===] made complete graph with {} edges".format(link_count))


# def make_infomap_initial(args, ctg_lookup: Dict[str, NbContig]) -> None:
#     def get_ctg2bin(bin_fa:str, ctg2bin_:Dict[int, int], bin_idx:int) -> None:
#         with open(bin_fa, 'r') as fr:
#             for iline in fr:
#                 iline = iline.strip()
#                 if iline.startswith(">"):
#                     iline = format_ctg_idf(iline.lstrip(">"))
#                     if iline in ctg_lookup:
#                         ctg2bin_[ctg_lookup[iline].num_id] = bin_idx
#     ctg2bin = {}
#     ibin_idx = 0
#     # need fix: make mapbin accept only one bin_dir
#
#     print("[=== Set initial clusters ===] Parsing binning result directory \"{}\"...".format(args.bin_dir))
#
#     bin_fa_files = glob.glob(args.bin_dir + "*." + args.bin_suffix)
#
#     for i_bin_fa in bin_fa_files:
#         ibin_idx += 1
#         get_ctg2bin(i_bin_fa, ctg2bin, ibin_idx)
#
#     args.init_clu = args.out_dir + "infomap.initial.clu"
#
#     print("[=== Set initial clusters ===] Write out binning-based initial clusters to \"{}\"".format(args.init_clu))
#     with open(args.init_clu, 'w') as fw:
#         line_maker = lambda k, v: "{} {}".format(k,v)
#         fw.write("\n".join([line_maker(k,v) for k, v in ctg2bin.items()]))


def _make_loners(loners: Set, im_network: Infomap, layer_idx:int, weight_factor) -> None:
    for l in loners:
        im_network.add_multilayer_intra_link(layer_idx, l, l, weight_factor)


# def _add_link_pairwise(x: NbContig, y: NbContig, bin_weighter:float,
#                        im_network:Infomap, bin_layer_idx:int):
#     # """long contig goes to short. short no go to long"""
#     """short contig goes to long. long no go to short"""
#     def _calc_edge_weight(tail: NbContig) -> float:
#         return sqrt(tail.length * bin_weighter)
#     if x.length < y.length:
#         w = _calc_edge_weight(y)
#         im_network.add_multilayer_intra_link(bin_layer_idx, x.num_id, y.num_id, w)
#     else:
#         w = _calc_edge_weight(x)
#         im_network.add_multilayer_intra_link(bin_layer_idx, x.num_id, y.num_id, w)


def make_infomap_binned(im_network: Infomap, layer: Dict[str, int],
                        ctg_lookup: Dict[str, NbContig],
                        bin_lookup: Dict[str, int]) -> None:
    layer_idx = layer['bin']

    def _calc_bin_weighter(n_bin:int) -> float:
        """
        scale the edges of the bin in the context of all bins
        :param n_bin: number of contigs in this bin
        :return: an efficient for the bin
        """
        # return 1/(sum_bases * n_bin * (n_bin - 1))
        # return sqrt(len_tail / (L_total * N_bin * (N_bin - 1)))
        return 1/total_bases*(2*n_bin - 2* sqrt(n_bin))

    def _calc_edge_weight(tail: NbContig, num_bases:int) -> float:
        return sqrt(tail.length /num_bases)

    def _add_link_pairwise(x: NbContig, y: NbContig, bin_weight:float):
    # def _add_link_pairwise(x: NbContig, y: NbContig):
        """long contig goes to short. short no go to long"""
        # if x.length > y.length:
        #     w = y.length * bin_weight  # short-len/total-contig-len*nlinks-in-this-bin
        #     # w = y.length/x.length
        #     # w = y.length/max_ctg_len
        #
        # else:
        #     w = x.length * bin_weight
        #     w = x.length/y.length
        #     w = x.length / max_ctg_len
        w = min(x.length, y.length)/(x.length + y.length)
        im_network.add_multilayer_intra_link(layer_idx, x.num_id, y.num_id, w)
        # im_network.add_multilayer_intra_link(layer_idx, x.num_id, y.num_id, w)
        # im_network.add_multilayer_intra_link(layer_idx, y.num_id, x.num_id, w)


    def _make_grid_row(ends_pool: List[NbContig], inner_pool: List[NbContig],
                       len_row: int, bin_weighter: float) -> List[NbContig]:
    # def _make_grid_row(ends_pool: Set[NbContig], inner_pool: Set[NbContig],
    #                    len_row: int) -> List[NbContig]:
        """
        General way to make a row of nodes, and add row edges.
        1st and last line: two ends are from corners, middle are from smalls.
        middle rows: two ends are from smalls, middle from rest.
        applies to both regular and incomplete grid graphs.
        """
        if len_row == 1:
            row = [ends_pool[0]]
            del ends_pool[0]
            return row

        elif len_row == 2:
            leftmost = ends_pool[0]
            rightmost = ends_pool[1]
            _add_link_pairwise(leftmost, rightmost, bin_weighter)
            # _add_link_pairwise(leftmost, rightmost)
            del ends_pool[:2]
            return [leftmost, rightmost]

        elif len_row > 2:
            row = []
            leftmost = ends_pool[0]
            rightmost = ends_pool[1]
            del ends_pool[:2]
            row.append(leftmost)
            for i in range(0, len_row - 2):
                inode = inner_pool[i]
                row.append(inode)  # inode's index in row is i + 1
                left_node = row[i]
                _add_link_pairwise(left_node, inode, bin_weighter)
                # _add_link_pairwise(left_node, inode)
            del inner_pool[:(len_row - 2)]

            row.append(rightmost)
            _add_link_pairwise(row[-2], rightmost, bin_weighter)
            # _add_link_pairwise(row[-2], rightmost)

            return row

    def _add_col_edges(grid: List[List[NbContig]], bin_weighter:float):
    # def _add_col_edges(grid: List[List[NbContig]]):

        for i in range(len(grid) - 1):  # nrow
            itop = grid[i]
            ibtm = grid[i+1]
            for j in range(len(ibtm)):
                x_node = itop[j]
                y_node = ibtm[j]
                _add_link_pairwise(x_node, y_node, bin_weighter)
            if len(ibtm) == 1:
                x_node = itop[1]
                y_node = ibtm[0]
                _add_link_pairwise(x_node, y_node, bin_weighter)
                # _add_link_pairwise(x_node, y_node)

    def make_bin_grid(sorted_nbctg: List[NbContig], bin_weighter:float):
    # def make_bin_grid(sorted_nbctg: List[NbContig]):
        N = len(sorted_nbctg)
        # if N == 1:
        #     im_network.addNode(sorted_nbctg[0].num_id, sorted_nbctg[0].id)
        if N == 2:
            _add_link_pairwise(sorted_nbctg[0], sorted_nbctg[1], bin_weighter)
            # _add_link_pairwise(sorted_nbctg[0], sorted_nbctg[1])
        elif 2 < N <= 5:
            for i in range(N - 1):
                x_node = sorted_nbctg[i]
                y_node = sorted_nbctg[i + 1]
                _add_link_pairwise(x_node, y_node, bin_weighter)
                # _add_link_pairwise(x_node, y_node)
            _add_link_pairwise(sorted_nbctg[-1], sorted_nbctg[0], bin_weighter)
            # _add_link_pairwise(sorted_nbctg[-1], sorted_nbctg[0])


        elif N > 5:
            ncol = ceil(sqrt(N))
            nrow = ceil(N/ncol)

            len_last = N - ncol * (nrow - 1)
            if len_last == 1:
                n_corners = 3
            else:
                n_corners = 4
            n_internal = (ncol - 2) * (nrow - 2)

            corners = sorted_nbctg[:n_corners]
            rest = sorted_nbctg[n_corners: N]

            # make the first row
            grid_net = [_make_grid_row(corners, rest, ncol, bin_weighter)]
            # grid_net = [_make_grid_row(corners, rest, ncol)]

            if nrow > 2:  # make middle rows
                for i in range(nrow - 2):
                    grid_net.append(_make_grid_row(rest, rest, ncol, bin_weighter))
                    # grid_net.append(_make_grid_row(rest, internal, ncol))

            line_last = _make_grid_row(corners, rest, len_last, bin_weighter)
            # line_last = _make_grid_row(corners, rest, len_last)
            grid_net.append(line_last)
            _add_col_edges(grid_net, bin_weighter)


            # corners = set(sorted_nbctg[:n_corners])
            # internal = set()
            # if n_internal:
            #     internal = set(sorted_nbctg[-(n_internal + 1): -1])
            # rest = set(sorted_nbctg[n_corners: (N - n_internal)])
            #
            # # make the first row
            # grid_net = [_make_grid_row(corners, rest, ncol, bin_weighter)]
            # # grid_net = [_make_grid_row(corners, rest, ncol)]
            #
            # if nrow > 2:
            #     for i in range(nrow - 2):
            #         grid_net.append(_make_grid_row(rest, internal, ncol, bin_weighter))
            #         # grid_net.append(_make_grid_row(rest, internal, ncol))
            #
            # line_last = _make_grid_row(corners, rest, len_last, bin_weighter)
            # # line_last = _make_grid_row(corners, rest, len_last)
            # grid_net.append(line_last)
            # _add_col_edges(grid_net, bin_weighter)
            # _add_col_edges(grid_net)

    # remaining_nodes = set([i.num_id for i in ctg_lookup.values()])
    binned_bins = set(bin_lookup.values())
    total_bases = sum([i.length for i in ctg_lookup.values()])

    for i in binned_bins:
        contigs = [ctg_lookup[ctg] for ctg, bin_id in bin_lookup.items() if bin_id == i]
        if len(contigs) == 1:
            continue
        print("[=== Binning layer ===] Make a grid-like network for {} contigs".format(len(contigs)))
        nbases = sum([i.length for i in contigs])
        # bw = _calc_bin_weighter(len(contigs))
        max_bp = max([i.length for i in contigs])
        bw = 0.5/max_bp

        contigs.sort(key=lambda x: x.length)
        make_bin_grid(contigs, bw)
            # make_bin_grid(contigs)


# def make_pairing_flat_net(args, ctg_lookup:Dict[str, NbContig]):
#     im_flat = Infomap("--undirected -s {}".format(args.im_seed))
#     time_begin = time.time()
#     print("[=== Info ===] Parsing alignment...")
#     ctg_links = process_sam(args.aln_file, args.aln_query_perc, ctg_lookup)
#     n_edges = 0
#     for (i, j), (aligned_ij, aligned_ji) in ctg_links.items():
#         if aligned_ij > args.min_ctg_len and aligned_ji > args.min_ctg_len:
#             ctg_i = ctg_lookup[i]
#             ctg_j = ctg_lookup[j]
#             im_flat.add_link(ctg_i.num_id, ctg_j.num_id,
#                              (aligned_ij + aligned_ji)/(ctg_i.length + ctg_j.length))
#             n_edges += 1
#
#     im_flat.run(ftree=True)
#     im_flat.write_pajek(args.im_pajek)
#     im_flat.write_flow_tree(args.im_ftree)
#     im_flat.write_clu(args.im_clu)
#
#     time_finish = time.time()
#     print("[=== 1st order network with read pairing ===] Made {} edges based on read pairing.".format(n_edges))
#     print("[=== 1st order network with read pairing ===] Processing read pairing took {:.1f} secs".format(
#         time_finish - time_begin))


def _add_mllink_from_ctg_links(ctg_i: NbContig, ctg_j:NbContig,
                               aligned_bp: int,
                               pair_layer_idx: int,
                               im_network: Infomap,
                               max_ctg_len:int):
    w_ij = aligned_bp/max_ctg_len
    # w_ji = aligned_ji / _w
    new_edges = 0
    # if w_ij > 0.001:
    im_network.add_multilayer_intra_link(pair_layer_idx, ctg_j.num_id,
                                         ctg_i.num_id, w_ij)
    new_edges += 1
    # if w_ji > 0.001:
    #     im_network.add_multilayer_intra_link(pair_layer_idx, ctg_j.num_id,
    #                                          ctg_i.num_id, w_ji)
    #     new_edges += 1
    return new_edges


def make_infomap_pair(args, im_multi: Infomap, layer: Dict[str, int],
                      ctg_lookup: Dict[str, NbContig],
                      bin_lookup: Dict[str, int],
                      n_process: int,
                      max_ctg_len: int):
    nproc = get_n_subprocesses(n_process, len(args.aln_files))
    ctg_links = process_bams(args.aln_files, args.aln_query_perc,
                             ctg_lookup, nproc, bin_lookup,
                             args.rlen, args.max_insert_size)
    pair_layer_idx = layer['pair']
    print("======= ctg_links:", len(ctg_links))

    for (ctg_i_id, ctg_j_id), aligned in ctg_links.items():
        ctg_i = ctg_lookup[ctg_i_id]
        ctg_j = ctg_lookup[ctg_j_id]
        ctg = sorted([ctg_i, ctg_j], key=lambda x: x.length)

        _add_mllink_from_ctg_links(ctg[0], ctg[1], aligned,
                                   pair_layer_idx, im_multi, max_ctg_len)




def make_infomap_assembly(args, im_network: Infomap,
                          layer: Dict[str, int],
                          ctg_lookup: Dict[str, NbContig],
                          bin_lookup: Dict[str, int],
                          max_ctg_len:int):
    layer_idx = layer['assembly']
    time_begin = time.time()
    seg_len_lookup = None
    links_lookup = None
    ctg2segs = None
    counter = 0
    # rejected_counter = 0
    # eff = log2(1/(max_ctg_len * max_ctg_len))

    if args.assembler == "spades":
        if not args.ctg_paths:
            seg_len_lookup, links_lookup, ctg2segs = parse_spades_gfa(ctg_lookup, args.gfa,
                                                                      parse_paths=True)
        else:
            seg_len_lookup, links_lookup, _ = parse_spades_gfa(ctg_lookup, args.gfa,
                                                               parse_paths=False)
            ctg2segs = parse_spades_contig_paths(ctg_lookup, args.ctg_paths)

        linked_ctg = link_spades_contigs(seg_len_lookup, links_lookup, ctg2segs,
                                         ctg_lookup, bin_lookup)
        for (i_ctg_name, j_ctg_name), linked_bp in linked_ctg.items():
            ic = ctg_lookup[i_ctg_name]
            jc = ctg_lookup[j_ctg_name]
            # ic_binned = -1
            # jc_binned = -1

            # if i_ctg_name in bin_lookup:
            #     ic_binned = bin_lookup[i_ctg_name]
            # if j_ctg_name in bin_lookup:
            #     jc_binned = bin_lookup[j_ctg_name]

            # if ic_binned == jc_binned or \
            #     ic_binned == -1 or \
            #         jc_binned == -1:
            #
            #     ctg = sorted([ic, jc], key=lambda x:x.length)
            #     ctg_lengths = [i.length if i.length < 10000 else 10000 for i in ctg ]
            #
            #     w = linked_bp*ctg[0].length/(ctg_lengths[0] + ctg_lengths[1] - 2*linked_bp)/max_ctg_len
            #     if w > 0.0005:
            #         print("==========assembly", (i_ctg_name, j_ctg_name), linked_bp, w )
            #         im_network.add_multilayer_intra_link(layer_idx, ctg[1].num_id,
            #                                      ctg[0].num_id, w)
            #
            #     counter += 1

            # ctg = sorted([ic, jc], key=lambda x: x.length)
            w = 2 * linked_bp / (ic.length + jc.length)

            # if ic_binned == jc_binned or \
            #     ic_binned == -1 or \
            #         jc_binned == -1:
                # ctg = sorted([ic, jc], key=lambda x:x.length)
                # w = 2* linked_bp/(ctg[0].length + ctg[1].length)
                # print(ctg[1].num_id, ctg[0].num_id, w)
                # w_ij = linked_bp/ic.length
                # w_ji = linked_bp/jc.length
                # print("==========assembly", (i_ctg_name, j_ctg_name), linked_bp )
                # if w > 0.01:
            print("==========assembly", (i_ctg_name, j_ctg_name), linked_bp, w )
            im_network.add_multilayer_intra_link(layer_idx, ic.num_id,
                                                 jc.num_id, w)
            im_network.add_multilayer_intra_link(layer_idx, jc.num_id,
                                                 ic.num_id, w)
            counter += 1
            # else:
            #     if w > 0.05:
            #         im_network.add_multilayer_intra_link(layer_idx, ic.num_id,
            #                                              jc.num_id, w)
            #         im_network.add_multilayer_intra_link(layer_idx, jc.num_id,
            #                                              ic.num_id, w)
            #     else:
            #         rejected_counter += 1

        print("[=== Assembly layer ===] made {} links from the assembly graph".format(counter))
        # print("[=== Info ===] rejected {} links because the two contigs are in different bins".format(rejected_counter))

    time_end = time.time()
    print("[=== Assembly layer ===] Processing assembly took {:.1f} secs".format(time_end - time_begin))


def generate_ctg_lookup(contig_fa_lookup: Dict[str, str], min_ctg_len: int):
    # alternative func for generate ctg lookup
    # function: generate ctg lookup dict
    #           and initialize im by initializing the network and adding nodes
    """
    :param contig_fa_lookup: {fasta_header: seq_len}
    :param min_ctg_len:
    :return: contig_fa_lookup.keys(): all contig names as a set
    :return: ctg_lookup: Dict[str, NbContig]
            contigs that passed the min contig length check, formatted for use in infomap
    :return: unbinned: [fasta record as string]
    :return: im_network is no longer generated here.
    """
    # im_network = Infomap("-s 1")
    ctg_lookup = {}
    unbinned_short = []
    i = 1

    for iheader, iseq in contig_fa_lookup.items():
        ilen = len(iseq)

        if ilen < min_ctg_len:
            unbinned_short.append(iheader)
        else:
            formatted_id = format_ctg_idf(iheader)
            cur_nbcontig = NbContig(formatted_id, ilen, i)  # changed node_id to int
            ctg_lookup[formatted_id] = cur_nbcontig  # changed the key to the formatted_id, mar 2022
            # im_network.set_name(i, formatted_id)
            # im_network.add_node(node_id=i,
            #                     name=formatted_id)

        i += 1
    print("[=== Parsing contigs ===] {} contigs with minimum length {} for binning".format(len(ctg_lookup), min_ctg_len))

    return ctg_lookup, unbinned_short


def init_infomap(args, ctg_lookup: Dict[str, NbContig]) -> Infomap:
    # if args.init_clu:
    #     im_network = Infomap("--variable-markov-time -s {}".format(args.im_seed), cluster_data=args.init_clu)
    # else:
    im_network = Infomap(
                         # "--prefer-modular-solution "
                         "--clu "
                         # "--variable-markov-damping 0.5 "
                         "--assign-to-neighbouring-module "
                         "--multilayer-relax-rate 0.3 "
                         # "--variable-markov-time "
                         # "--core-loop-codelength-threshold 0.001 "
                         # "--tune-iteration-relative-threshold 0.1 "
                         "--teleportation-probability 0.15 ", num_trials=10)
                         # "--teleportation-probability 0.3 -s {} ".format(args.im_seed), num_trials=100)

    for k, nbctg in ctg_lookup.items():
        im_network.set_name(nbctg.num_id, k)

    return im_network
