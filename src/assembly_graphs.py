from typing import Set, List, Tuple, Dict, Iterable
from seq_objects import *
import sys
from statistics import median
from collections import Counter, namedtuple


def parse_spades_contig_paths(ctg_lookup: Dict[str, NbContig],
                              path_file: str) -> Dict[int, Tuple[Set[str], Set[str]]]:
    """
    :param: ctg_lookup: Dict[str, NbContig]. contains contigs > min_ctg_len.
            keys are full fasta header
            values are formatted info to be used in infomap
    :param: all_contigs: keys from contig_fa_lookup: Dict[str, str].
    :returns: {contig: [(path_left_end, path_right_end)]}, contains only the contigs whose length > min_ctg_len
    :returns: {contig_names: (<list of path segments>, <list of ends>)}
    """

    def _update_ctg2segs(ctg_name: str, path_segs: List[str]):
        # path_segs contained list of RAW path strings (semicomma removed
        if path_segs:
            for p in path_segs:
                ipath = [i[:-1] for i in p.split(",")]
                ends = {ipath[0], ipath[-1]}
                ctg_numid = ctg_lookup[ctg_name].num_id
                if ctg_name in ctg2segs:
                    ctg2segs[ctg_numid][0].update(ipath)
                    ctg2segs[ctg_numid][1].update(ends)
                else:
                    ctg2segs[ctg_numid] = (set(ipath), ends)

    ctg2segs = {}
    with open(path_file, 'r') as fr:
        is_rev = 0
        cur_header = ""
        cur_rev_header = ""
        cur_path = []
        path_found = False
        for iline in fr:
            iline = iline.rstrip()
            if iline:
                if iline in ctg_lookup:
                    if not path_found:
                        path_found = True
                    cur_header = iline
                    is_rev = 0
                    cur_rev_header = iline + "\'"

                elif iline == cur_rev_header:
                    is_rev = 1
                    # conclude the current entry
                    if cur_path:
                        # if cur_header in ctg_lookup:
                        _update_ctg2segs(cur_header, cur_path)
                    cur_path = []
                    # cur_header = ""

                elif is_rev:
                    continue
                else:
                    cur_path.append(iline.rstrip(";"))

    print("[=== Assembly layer ===] Parsed assembly paths from path file: {}".format(path_file))
    if not path_found:
        print("[=== Assembly layer ===] No paths found to match contig in "
              "the fasta file. Sure the files match?")

    return ctg2segs


def parse_spades_gfa(ctg_lookup: Dict[str, NbContig],
                     gfa_file: str, parse_paths: bool):
    def _parse_lcigar(cig: str):
        n_match = 0
        try:
            n_match = int(cig[:-1])
        except ValueError:
            pass
        return n_match

    def _parse_s_line(sline: str):
        sline_split = sline.split("\t")
        _, sid, segseq = sline_split[:3]

        return sid, len(segseq), int(sline_split[-1].split(":")[-1])

    def _update_ctg2segs(ctg_name: str, path_segs: List[str]):
        ends = {path_segs[0], path_segs[-1]}
        ctg_numid = ctg_lookup[ctg_name].num_id
        if ctg_numid in ctg2segs:
            ctg2segs[ctg_numid][0].update(path_segs)
            ctg2segs[ctg_numid][1].update(ends)
        else:
            ctg2segs[ctg_numid] = (set(path_segs), ends)

    ctg2segs = {}  # str: Tuple[Set[str],Set[str]]]
    seg_lookup = {}  # str(seg-name): List[seg_len, depth]
    links_lookup = {}  # tuple[str, str]: int(overlap-len)
    path_found = False
    with open(gfa_file, 'r') as fr:
        for iline in fr:
            iline = iline.rstrip()
            if iline:
                if iline.startswith("S"):
                    seg_name, seg_len, kc = _parse_s_line(iline)
                    seg_lookup[seg_name] = (seg_len, kc / seg_len, [])

                elif iline.startswith("L"):
                    _, i, _, j, _, lap = iline.split("\t")[:6]
                    nlap = _parse_lcigar(lap)

                    if nlap:
                        links_lookup[(i, j)] = nlap

                elif iline.startswith("P"):
                    if not parse_paths:
                        continue
                    _, path_id, segments = iline.split("\t")[:3]
                    ictg = "_".join(path_id.split("_")[:-1])
                    if ictg in ctg_lookup:
                        segments = [i[:-1] for i in segments.split(",")]
                        _update_ctg2segs(ictg, segments)
                        if not path_found:
                            path_found = True
    print("[=== Assembly layer ===] Parsed assembly graph from gfa file: {}".format(gfa_file))
    if parse_paths and not path_found:
        print("[=== Assembly layer ===] No paths found to match contig in "
              "the fasta file. Sure the files match?")

    return seg_lookup, links_lookup, ctg2segs


def _calc_ctg_depth(ctg2segs: Dict[int, Tuple[Set[str], Set[str]]],
                    seg_lookup: Dict[str, Tuple[int, float, List[int]]]):
    ctg_depth = {}
    # generate ctg_depth lookup
    try:
        for k, v in ctg2segs.items():
            ctg_depth[k] = median([seg_lookup[iseg][1] for iseg in v[0]])
    except KeyError:
        print("[=== Error ===] Segment not found in graph. The graph file and path file may not match",
              file=sys.stderr)
        sys.exit(1)

    return ctg_depth


def _add_ctg_names_2_seglookup(seg_lookup: Dict[str, Tuple[int, float, List[int]]],
                               ctg2segs: Dict[int, Tuple[Set[str], Set[str]]]):
    """it appends the contigs to seglookup"""
    for numid, (segs, _) in ctg2segs.items():
        for _iseg in segs:
            seg_lookup[_iseg][2].append(numid)


def _link_spades_contigs_shared_segs(seg_lookup: Dict[str, Tuple[int, float, List[int]]],
                                     ctg_depth: Dict[int, float]):
    """ Link contigs by their shared segments. also remove
    contigs whose depth is way below that of segment """

    # update shared_bp
    ctg_links = {}
    new_seg_lookup = {}
    for seg_name, (s_len, s_dp, s_ctg) in seg_lookup.items():
        if len(s_ctg) < 2:
            new_seg_lookup[seg_name] = (s_dp, seg_lookup[seg_name][2])
            continue

        kept_ctg = [i for i in s_ctg if s_dp < 10 * ctg_depth[i]]
        len_kept_ctg = len(kept_ctg)
        new_seg_lookup[seg_name] = (s_dp, kept_ctg)

        if len_kept_ctg < 2:
            continue

        kept_ctg.sort()
        for i in range(len_kept_ctg):
            for j in range(i + 1, len_kept_ctg):
                ictg_numid = kept_ctg[i]
                jctg_numid = kept_ctg[j]
                if (ictg_numid, jctg_numid) in ctg_links:
                    ctg_links[(ictg_numid, jctg_numid)] += s_len
                else:
                    ctg_links[(ictg_numid, jctg_numid)] = s_len

    return new_seg_lookup, ctg_links


def link_spades_contigs(seg_lookup: Dict[str, Tuple[int, float, List[int]]],
                        links_lookup: Dict[Tuple[str, str], int],
                        ctg2segs: Dict[int, Tuple[Set[str], Set[str]]],
                        ctg_lookup: Dict[str, NbContig],
                        bin_lookup: Dict[str, int]) -> Dict[Tuple[str, str], int]:
    linkmate = namedtuple("LinkMate", ["node_id", "node_bin", "linked_bp"])

    def _get_bin_id(ctg_id: str):
        if ctg_id in bin_lookup:
            return bin_lookup[ctg_id]
        else:
            return -1

    def _update_ctg_links(x: int, y: int, bp: int):
        key = tuple(sorted([x, y]))
        if key in _ctg_links:
            _ctg_links[key] += bp
        else:
            _ctg_links[key] = bp

    def filter_links(links: Dict[Tuple[int, int], int]):
        """filter links and change the keys from int to str names"""

        def _update_nodes2nodes(x: str, x_linkmate):
            if x in _nodes2nodes:
                _nodes2nodes[x].append(x_linkmate)
            else:
                _nodes2nodes[x] = [x_linkmate]


        _to_remove = set()
        _nodes2nodes = {}
        too_small_count = 0
        ctg_links = {}

        for (k0, k1), v in links.items():
            k0_id = _ctg_numid_lookup[k0]
            k1_id = _ctg_numid_lookup[k1]
            if v < 300:
                too_small_count += 1
                continue

            k0_bin = _get_bin_id(k0_id)
            k1_bin = _get_bin_id(k1_id)

            if k0_bin == -1 and v < 0.1 * ctg_lookup[k0_id].length:
                continue
            elif k1_bin == -1 and v < 0.1 * ctg_lookup[k1_id].length:
                continue

            k0_linkmate = linkmate(k1_id, k1_bin, v)
            k1_linkmate = linkmate(k0_id, k0_bin, v)
            _update_nodes2nodes(k0_id, k0_linkmate)
            _update_nodes2nodes(k1_id, k1_linkmate)

            ctg_links[(k0_id, k1_id)] = v

        for i, i_linked in _nodes2nodes.items():
            i_bin = _get_bin_id(i)

            if len(i_linked) == 1:
                il_bin = i_linked[0].node_bin
                if i_bin == -1 or il_bin == -1 or i_bin == il_bin:
                    continue
                else:
                    if i_linked[0].linked_bp > ctg_lookup[i].length:
                        print(">>>>>>>>>>>>> removing from bin {}, {}".format(i_bin, i))
                        del bin_lookup[i]
                    else:
                        _to_remove.update({(i, i_linked[0].node_id), (i_linked[0].node_id, i)})
            else:
                bin_counter = Counter()
                for il in i_linked:
                    if il.node_bin != -1:
                        if il.node_bin in bin_counter:
                            bin_counter[il.node_bin] += il.linked_bp
                        else:
                            bin_counter[il.node_bin] = il.linked_bp

                kept_links = []
                if len(bin_counter) != 0:
                    largest_bin = bin_counter.most_common(1)[0]
                    kept_bin = largest_bin[0]
                    if i_bin != -1:
                        # if the binned contig is more linked to other bins:
                        if bin_counter[i_bin] == 0 and largest_bin[1] > ctg_lookup[i].length:
                            print(">>>>>>>>>>>>> removing from bin {}, {}".format(i_bin, i))
                            del bin_lookup[i]
                        else:
                            kept_bin = i_bin

                    for il in i_linked:
                        if il.node_bin != kept_bin and il.node_bin != -1:
                            _to_remove.update({(i, il.node_id), (il.node_id, i)})
                        else:
                            kept_links.append(il)

                else:
                    kept_links = i_linked

                if len(kept_links) > 1:
                    i_max_bp = max([k.linked_bp for k in kept_links])
                    for il in kept_links:
                        if il.linked_bp < 0.5 * i_max_bp:
                            _to_remove.add((i, il.node_id))

        for (k0, k1) in _to_remove:
            if (k1, k0) not in _to_remove:
                continue
            if (k0, k1) in ctg_links:
                del ctg_links[(k0, k1)]
            elif (k1, k0) in ctg_links:
                del ctg_links[(k1, k0)]
        print("=======too small", too_small_count)

        return ctg_links

    _add_ctg_names_2_seglookup(seg_lookup, ctg2segs)
    ctg_depth = _calc_ctg_depth(ctg2segs, seg_lookup)
    new_seg_lookup, _ctg_links = _link_spades_contigs_shared_segs(seg_lookup, ctg_depth)

    for (i_seg, j_seg), lapped_bp in links_lookup.items():
        # contigs sharing the segment can be linked
        i_depth, ictgs = new_seg_lookup[i_seg]
        j_depth, jctgs = new_seg_lookup[j_seg]
        if ictgs and jctgs:
            for iictg in ictgs:
                for jjctg in jctgs:
                    if iictg != jjctg:
                        _depths = sorted([ctg_depth[iictg], ctg_depth[jjctg], i_depth, j_depth])
                        if max(i_depth, j_depth) < 10 * min(ctg_depth[iictg], ctg_depth[jjctg]):
                            _update_ctg_links(iictg, jjctg, lapped_bp)

    _ctg_numid_lookup = {v.num_id: k for k, v in ctg_lookup.items()}

    ctg_links = filter_links(_ctg_links)

    # _nodes2nodes = {}  # node_ids to a list of node_ids that they link to
    # _to_remove = set()
    # _ctg_links = {}  # node_ids to linked_bp
    # for (k0, k1), v in ctg_links.items():
    #     # filtered short links, between bin short links
    #     if v < 300:
    #         continue
    #     k0_id = _ctg_numid_lookup[k0]
    #     k1_id = _ctg_numid_lookup[k1]
    #
    #     k0_bin = _get_bin_id(k0_id)
    #     k1_bin = _get_bin_id(k1_id)
    #
    #     if k0_bin != -1 and k1_bin != -1 \
    #             and k0_bin != k1_bin:
    #         k0_perc = v/ctg_lookup[k0_id].length
    #         k1_perc = v/ctg_lookup[k1_id].length
    #         _taken = False
    #         if k0_perc > 1/3 or (v > 2000 and k0_perc > 0.05):
        #         k0_bin = -1
        #         print(k0_id, k1_id)
        #         del bin_lookup[k0_id]
        #         _taken = True
        #     if k1_perc > 1/3 or (v > 2000 and k1_perc > 0.05):
        #         k1_bin = -1
        #         print(k0_id, k1_id)
        #         del bin_lookup[k1_id]
        #         _taken = True
        #
        #     # _smaller = min(ctg_lookup[k0_id], ctg_lookup[k1_id], key=lambda x:x.length)
        #     # if v/_smaller.length > 1/3 or \
        #     #         (v > 2000 and v/_smaller.length > 0.05):
        #     #     print(k0_id, k1_id)
        #     #     bin_lookup.pop(_smaller.id)
        #     if _taken:
        #         k0_linkmate = linkmate(k1_id, k1_bin, v)
        #         _update_nodes2nodes(k0_id, k0_linkmate)
        #         k1_linkmate = linkmate(k0_id, k0_bin, v)
        #         _update_nodes2nodes(k1_id, k1_linkmate)
        #         _ctg_links[tuple(sorted([k0_id, k1_id]))] = v
        #
        # else:
        #     k0_linkmate = linkmate(k1_id, k1_bin, v)
        #     _update_nodes2nodes(k0_id, k0_linkmate)
        #     k1_linkmate = linkmate(k0_id, k0_bin, v)
        #     _update_nodes2nodes(k1_id, k1_linkmate)
        #     _ctg_links[tuple(sorted([k0_id, k1_id]))] = v

    # for k, k_linked in _nodes2nodes.items():
    #     if len(k_linked) == 1:
    #         # print("==== single link", k, k_linked[0].node_id)
    #         # _ctg_links[tuple(sorted([k, k_linked[0].node_id]))] = k_linked[0].linked_bp
    #         continue
    #     k_bin = _get_bin_id(k)
    #     k_max_bp = max([i.linked_bp for i in k_linked])
    #     # total_bp = 0
    #     k_kept = []
    #     for i in k_linked:
    #         if i.linked_bp < 0.5 * k_max_bp:
    #             _to_remove.add(tuple(sorted([k, i.node_id])))
    #         else:
    #             k_kept.append(i)
    #
    #     if k_bin == -1:
    #         bin_counter = Counter()
    #         for i in k_kept:
    #             # total_bp += i.linked_bp
    #             if i.node_bin != -1:
    #                 if i.node_bin in bin_counter:
    #                     bin_counter[i.node_bin] += i.linked_bp
    #                 else:
    #                     bin_counter[i.node_bin] = i.linked_bp
    #         if len(bin_counter) < 2:
    #             # print("===== One bin")
    #             # for kl in k_linked:
    #             #     print(k, kl.node_id)
    #             #     _ctg_links[tuple(sorted([k, kl.node_id]))] = k_linked[0].linked_bp
    #             continue
    #
    #         largest_bin = bin_counter.most_common(1)[0]
    #         if largest_bin[1] >= k_max_bp:
    #             print("===== multiple bin removed")
    #             for kl in k_kept:
    #                 if kl.node_bin != largest_bin[0] and kl.node_bin != -1:
    #                     print(k, kl.node_id, kl.linked_bp)
    #                     _to_remove.add(tuple(sorted([k, kl.node_id])))
    #
    #                     # _ctg_links[tuple(sorted([k, kl.node_id]))] = k_linked[0].linked_bp
    #         else:
    #             print("==== bins all removed")
    #             for kl in k_kept:
    #                 if kl.node_bin != -1:
    #                     _to_remove.add(tuple(sorted([k, kl.node_id])))
    #                     print(k, kl.node_id, kl.linked_bp)
    #                     # _ctg_links[tuple(sorted([k, kl.node_id]))] = k_linked[0].linked_bp
    #
    #     # else:
    #     #     print("==== binned ")
    #     #     for kl in k_linked:
    #     #         print(k, kl.node_id)
    #     #         _ctg_links[tuple(sorted([k, kl.node_id]))] = k_linked[0].linked_bp
    # for k in _to_remove:
    #     # print(k, _ctg_links[k])
    #     del _ctg_links[k]

    return ctg_links



    #         k0_linked = v / ctg_lookup[k0_id].length
    #         k1_linked = v / ctg_lookup[k1_id].length
    #
    #         if k0_linked > 0.1:
    #             k0_linkmate = linkmate(k1_id, k1_bin_id, v)
    #             _update_nodes2nodes(k0_id, k0_linkmate)
    #         if k1_linked > 0.1:
    #             k1_linkmate = linkmate(k0_id, k0_bin_id, v)
    #             _update_nodes2nodes(k1_id, k1_linkmate)
    #
    #     else:
    #         k0_linkmate = linkmate(k1_id, k1_bin_id, v)
    #         k1_linkmate = linkmate(k0_id, k0_bin_id, v)
    #         _update_nodes2nodes(k0_id, k0_linkmate)
    #         _update_nodes2nodes(k1_id, k1_linkmate)
    #
    # for k, k_linked in _nodes2nodes.items():
    #     k_bin = _get_bin_id(k)
    #     _kept = []
    #     _binned = []
    #     k_max_bp = max([i.linked_bp for i in k_linked])
    #     # keep at most one bin link for the node
    #
    #     k_linked_bins = []
    #     for i in k_linked:
    #         if i.linked_bp < 0.5 * k_max_bp:
    #             continue
    #
    #         if i.node_bin != -1:
    #             _binned.append(i)
    #             k_linked_bins.append(i.node_bin)
    #         else:
    #             _kept.append(i.node_id)
    #
    #     most_linked_bin = max(k_linked_bins, key=k_linked_bins.count)
    #     if k_bin != -1:
    #         if most_linked_bin != k_bin:
    #             del bin_lookup[k]
    #
    #     _kept.extend([i.node_id for i in _binned if i.node_bin == most_linked_bin])
    #
    # #
    # #
    # #
    # #     bin_max_bp = max(bin2bp.values())
    # #     if k_bin != -1 and bin2bp[k_bin] != bin_max_bp:
    # #         del bin_lookup[k]
    # #
    # #     elif
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #     k_max_bp = max([i.linked_bp for i in k_linked])
    # #
    # #     if k_bin != -1:
    # #         for i in k_linked:
    # #
    # #
    # #
    # #
    # #     if k_bin != -1:  # check if k should be unbinned
    # #         linked_same_bin_bp = sum([i.linked_bp for i in k_linked if i.bin_id == k_bin])
    # #         if linked_same_bin_bp < 0.5 * k_max_bp:
    # #             del bin_lookup[k]
    # #
    # #         for i in k_linked:
    # #             if i.bin_id == k_bin:
    # #                 linked_same_bin_bp += i.linked_bp
    # #             elif i.bin_id !=
    # #
    # #
    # #
    # #
    # #         linked_same_bin_bp = 0
    # #         for i in k_linked:
    # #             if i.bin_id == k_bin:
    # #                 linked_same_bin_bp += i.linked_bp
    # #         sum([i. for i in k_linked if i.bin_id == k_bin])
    # #
    # #
    # #         for i in range(nlinked):
    # #             if k_bin != k_linked[i].node_bin:
    # #                 iperc = k_linked[i].linked_bp/k_len
    # #
    # #     for i in range(nlinked):
    # #
    # #
    # #     if k_bin != -1:
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # # for (k0, k1), v in ctg_links.items():
    # #     if v < 300:
    # #         continue
    # #     k0_id = _ctg_numid_lookup[k0]
    # #     k1_id = _ctg_numid_lookup[k1]
    # #     key = tuple(sorted([k0_id, k1_id]))
    # #     if k0_id in bin_lookup and k1_id in bin_lookup:
    # #         k0_bin = bin_lookup[k0_id]
    # #         k1_bin = bin_lookup[k1_id]
    # #         if k0_bin != k1_bin:
    # #             k0_perc = v / ctg_lookup[k0_id].length
    # #             k1_perc = v / ctg_lookup[k1_id].length
    # #             if k0_perc > 0.1:
    # #                 del bin_lookup[k0_id]
    # #                 _ctg_links[key] = v
    # #             elif k1_perc > 0.1:
    # #                 del bin_lookup[k1_id]
    # #                 _ctg_links[key] = v
    # #         else:
    # #             _ctg_links[key] = v
    # #     else:
    # #         _ctg_links[key] = v
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #
    # #     if v >= 300:
    # #         _ctg_links[key] = v
    # #         _update_nodes2nodes(k0_id, k1_id)
    # #
    # # for k, k_linked in _nodes2nodes.items():
    # #     nlinked = len(k_linked)
    # #     if nlinked == 1:
    # #         kl = k_linked[0]
    # #         k_kl = tuple(sorted([k, kl]))
    # #         linked_bp = _ctg_links[k_kl]
    # #         if k in bin_lookup and kl in bin_lookup:
    # #             if bin_lookup[k] != bin_lookup[kl]:
    # #                 perc_k = linked_bp/ctg_lookup[k].length
    # #                 perc_kl = linked_bp/ctg_lookup[kl].length
    # #                 if max(perc_kl, perc_k) < 0.1:
    # #                     del _ctg_links[k_kl]
    # #                 elif perc_k >= 0.1:
    # #                     del bin_lookup[k]
    # #
    # #
    # #         # kl = k_linked[0]
    # #         # if k in bin_lookup and kl in bin_lookup:
    # #         #     if bin_lookup[k] != bin_lookup[kl]:
    # #         #         _to_remove.add(tuple(sorted([k, kl])))
    # #         continue
    # #     else:
    # #         linked_bins = set()
    # #         if k in bin_lookup:
    # #             linked_bins.add(bin_lookup[k])
    # #         else:
    # #             linked_bins.add(-1)
    # #
    # #         links = [tuple(sorted([k, kl])) for kl in k_linked]
    # #         weights = [_ctg_links[i] for i in links]
    # #         max_bp = max(weights)
    # #         _nlinked_filtered = 0
    # #         for i in range(nlinked):
    # #             if weights[i] < max_bp / 2:
    # #                 _to_remove.add(links[i])
    # #             else:
    # #                 _nlinked_filtered += 1
    # #                 kl = k_linked[i]
    # #                 if kl in bin_lookup:
    # #                     linked_bins.add(bin_lookup[kl])
    # #
    # #         if len(linked_bins) > 2 or _nlinked_filtered > 5:
    # #
    # #             _to_remove.update(links)
    # #         # else:
    # #             # print(k, weights, _nlinked_filtered, linked_bins)
    # #
    # #         # links = [links[i] for i in range(nlinked) if weights[i] > max_bp/3]
    # #
    # #     # linked_bins = set()
    # #     # if k in bin_lookup:
    # #     #     linked_bins.add(bin_lookup[k])
    # #     # for kl in k_linked:
    # #     #     if kl in bin_lookup:
    # #     #
    # #     #         linked_bins.add(bin_lookup[kl])
    # #
    # #
    # #
    # #
    # # for k in _to_remove:
    # #     # print(k, _ctg_links[k])
    # #     del _ctg_links[k]
    #
    # # _multiplicity_counter = Counter([ic for p in ctg_links.keys() for ic in p])
    #
    # # for (k0, k1), v in _ctg_links.items():
    # #     # if _multiplicity_counter[k0] > 3 or _multiplicity_counter[k1] > 3:
    # #     #     continue
    # #     # key = (_ctg_numid_lookup[k0],_ctg_numid_lookup[k1])
    # #     # _n = sorted([ctg_lookup[key[0]].length, ctg_lookup[key[1]].length])
    # #     _n = sorted([ctg_lookup[k0].length, ctg_lookup[k1].length])
    # #
    # #     if v > _n[0]:
    # #         v = _n[0]
    # #     # if v > l_perc * _n[0] and v > 0.1 * l_perc*_n[1]:
    # #     # if v> max(300, l_perc * _n[0], 0.1 * l_perc*_n[1]):
    # #     # if v > 0.005 * _n[1]:
    # #
    # #     new_ctg_links[key] = v


