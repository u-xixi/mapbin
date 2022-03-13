from infomap import Infomap
from typing import Dict, List, Tuple, Set


def fasta_concat(header: str, seq:str, line_len=70):
    """
    fasta file handler. format sequence line as fixed-length.
    """
    outstr = [">{}".format(header)]
    for i in range(len(seq)//line_len):
        outstr.append(seq[i*line_len: (i+1)*line_len])
    outstr.append(seq[(len(seq)//line_len)*line_len:])
    return "\n".join(outstr)


def _write_bin_dict(bin_dict:Dict[int, List[str]], out_dir: str):
    bin_dir = out_dir + "bins/"
    print("[=== Contig binning ===] writing binning output to {}".format(bin_dir))
    for ibin, icontigs in bin_dict.items():
        ibin_fa_filepath = "{}bin_{}.fa".format(bin_dir, ibin)
        with open(ibin_fa_filepath, 'w') as fw_handle:
            fw_handle.write("\n".join(icontigs) + "\n")


def _write_infomap_unbinned_contigs(binned_nodes: Set[str], contig_fa_lookup: Dict[str, str], out_dir:str) -> None:
    unbinned_content = []
    remains = set(contig_fa_lookup.keys()).difference(binned_nodes)

    for i in remains:
        unbinned_content.append(fasta_concat(header=i, seq=contig_fa_lookup[i]))
    unbinned_fp = out_dir + "bins/unbinned.fa"
    print("[=== Contig binning ===] {} contigs not assigned".format(len(unbinned_content)))

    if unbinned_content:
        print("[=== Contig binning ===] writing unbinned contigs to {}".format("unbinned.fa"))
        with open(unbinned_fp, 'w') as fw_handle:
            fw_handle.write("\n".join(unbinned_content) + "\n")


def write_unbinned_short(unbinned: List[str], contig_fa_lookup: Dict[str, str], outfile:str) -> None:
    """
    for unbinned contigs.
    should be used before writing out the infomap binned.
    """
    unbinned_short_content = []
    for iheader in unbinned:
        iseq = contig_fa_lookup.pop(iheader)
        unbinned_short_content.append(fasta_concat(iheader, iseq))

    with open(outfile, 'w') as fw:
        fw.write("\n".join(unbinned_short_content) + "\n")


def net2bin_overlap(im_network:Infomap, contig_fa_lookup: Dict[str, str], out_dir:str) -> None:
    """
    write out bins. bins might share contigs.
    :param contig_fa_lookup: contig_fa_lookup but with the contigs whose length < min_ctg_len have been popped out
    :param out_dir: args.out_dir  + bins/
    :return: None.
    write out:
    bin_x.fa files
    unbinned.fa <unbinned by infomap algo>
    disputed.tsv
    """
    def _pick_the_overlapped() -> Tuple[Dict[str, List[int]], Dict[str, List[int]]]:
        """
        all nodes are in nodes2modules, with all their modules
        if they are assigned to more than 1 module,
        they will also be in disputed.
        """
        nodes2modules = {}  # contig name: (module id, flow value)
        disputed = {}
        for i_node in im_network.nodes:
            i_module_id = i_node.module_id
            i_name = im_network.get_name(i_node.node_id)
            if i_name in nodes2modules:
                if i_name in disputed:
                    disputed[i_name].append(i_module_id)
                else:
                    disputed[i_name] = [nodes2modules[i_name][0], i_module_id]
                nodes2modules[i_name].append(i_module_id)

            else:
                nodes2modules[i_name] = [i_module_id]
        

        print("[=== Contig binning ===] {} contigs in more than 1 bins (disputed contigs)".format(len(disputed)))
        return nodes2modules, disputed

    def _write_disputed(disputed: Dict[str, List[int]]) -> None:
        if disputed:
            get_tabular_line = lambda k,v: "{}\t{}".format(k, ",".join([str(i) for i in v]))
            tabular_content = [get_tabular_line(k,v) for k, v in disputed.items()]
            unbinned_disputed_txt_fp = out_dir + "disputed.tsv"
            print(unbinned_disputed_txt_fp)

            print("[=== Contig binning ===] Report disputed contigs in {}".format("disputed.tsv"))
            with open(unbinned_disputed_txt_fp, 'w') as fw:
                fw.write("\n".join(tabular_content) + "\n")

    print("[=== Contig binning ===] binned into {} clusters".format(im_network.num_top_modules))
    print("[=== Contig binning ===] num_non_trivial_top_modules {}".format(im_network.num_non_trivial_top_modules))
    print("[=== Contig binning ===] num_effective_num_top_modules {}".format(im_network.effective_num_top_modules))
    print("[=== Contig binning ===] No. states: {}".format(im_network.num_nodes))
    print("[=== Contig binning ===] No. physical nodes: {}".format(im_network.num_physical_nodes))

    all_nodes_, disputed_ = _pick_the_overlapped()
    fa_content = {}  # module_id: fasta

    for header, bins in all_nodes_.items():
        for ib in bins:
            seq = contig_fa_lookup[header]
            if ib in fa_content:
                fa_content[ib].append(fasta_concat(header, seq))
            else:
                fa_content[ib] = [fasta_concat(header, seq)]

    _write_bin_dict(fa_content, out_dir)
    _write_infomap_unbinned_contigs(set(all_nodes_.keys()), contig_fa_lookup, out_dir)
    _write_disputed(disputed_)


def net2bin_nodispute(im_network:Infomap, contig_fa_lookup: Dict[str, str], out_dir:str) -> None:
    def _pick_the_overlapped() -> Tuple[Dict[str, int], Dict[str, List[int]]]:
        """
        Nodes are in nodes2modules, if they are assigned to just 1 module.
        otherwise they will be in disputed.
        """
        nodes2modules = {}
        disputed = {}
        for i_node in im_network.nodes:
            i_module_id = i_node.module_id
            i_name = im_network.get_name(i_node.node_id)
            if i_name in nodes2modules and i_name not in disputed:
                disputed[i_name] = [nodes2modules.pop(i_name)]
                disputed[i_name].append(i_module_id)
            elif i_name in disputed:
                disputed[i_name].append(i_module_id)
            else:
                nodes2modules[i_name] = i_module_id
        print("[=== Contig binning ===] {} contigs end up in more than 1 bins (disputed contigs)".format(len(disputed)))
        return nodes2modules, disputed

    def _write_disputed(disputed: Dict[str, List[int]]):
        if disputed:
            get_tabular_line = lambda k,v: "{}\t{}".format(k, ",".join([str(i) for i in v]))
            fasta_content = []
            tabular_content = []
            for k, v in disputed.items():
                tabular_content.append(get_tabular_line(k,v))
                k_seq = contig_fa_lookup[k]
                fasta_content.append(fasta_concat(k, k_seq))

            unbinned_disputed_txt_fp = out_dir + "disputed.tsv"
            unbinned_disputed_fa_fp = out_dir + "bins/unbinned.disputed.fa"

            print("[=== Contig binning ===] Report disputed contigs in {}".format("disputed.tsv"))
            with open(unbinned_disputed_txt_fp, 'w') as fw:
                fw.write("\n".join(tabular_content) + "\n")

            print("[=== Contig binning ===] Write disputed contigs as fasta to {}".format("unbinned.disputed.fa"))
            with open(unbinned_disputed_fa_fp, 'w') as fw:
                fw.write("\n".join(fasta_content) + "\n")

    print("[=== Contig binning ===] binned into {} clusters".format(im_network.num_top_modules))
    print("[=== Contig binning ===] No. states: {}".format(im_network.num_nodes))
    print("[=== Contig binning ===] No. physical nodes: {}".format(im_network.num_physical_nodes))

    undisputed_, disputed_ = _pick_the_overlapped()
    print(disputed_)
    undisputed_content = {}
    for i_header, i_bin in undisputed_.items():
        i_seq = contig_fa_lookup[i_header]
        if i_bin in undisputed_content:
            undisputed_content[i_bin].append(fasta_concat(i_header, i_seq))
        else:
            undisputed_content[i_bin] = [fasta_concat(i_header, i_seq)]


    _write_bin_dict(undisputed_content, out_dir)
    unbinned_infomap = set(contig_fa_lookup.keys()).difference(undisputed_.keys()).difference(disputed_.keys())
    _write_infomap_unbinned_contigs(unbinned_infomap, contig_fa_lookup, out_dir)
    _write_disputed(disputed_)


def net2bin_basic(im_network:Infomap, contig_fa_lookup: Dict[str, str], out_dir:str):
    def _pick_the_overlapped() -> Tuple[Dict[str, Tuple[int, float]], Dict[str, List[int]]]:
        """
        all nodes are in nodes2modules.
        if they are assigned to more than 1 module,
        they will also be in disputed.
        Their module in nodes2modules is the one that has the max flow
        """
        nodes2modules = {}  # contig name: (module id, flow value)
        disputed = {}
        for i_node in im_network.nodes:
            i_module_id = i_node.module_id
            i_name = im_network.get_name(i_node.node_id)
            i_flow = i_node.flow
            if i_name in nodes2modules:
                i_prev = nodes2modules[i_name]
                if i_flow > i_prev[1]:
                    nodes2modules[i_name] = (i_module_id, i_flow)
                if i_name in disputed:
                    disputed[i_name].append(i_module_id)
                else:
                    disputed[i_name] = [i_prev[0], i_module_id]

            else:
                nodes2modules[i_name] = (i_module_id, i_flow)

        print("[=== Contig binning ===] {} contigs in more than 1 bins (disputed contigs)".format(len(disputed)))
        return nodes2modules, disputed

    def _write_disputed(disputed: Dict[str, List[int]]) -> None:
        if disputed:
            get_tabular_line = lambda k,v: "{}\t{}".format(k, ",".join([str(i) for i in v]))
            tabular_content = [get_tabular_line(k,v) for k, v in disputed.items()]
            unbinned_disputed_txt_fp = out_dir + "disputed.tsv"

            print("[=== Contig binning ===] Report disputed contigs in {}".format("disputed.tsv"))
            with open(unbinned_disputed_txt_fp, 'w') as fw:
                fw.write("\n".join(tabular_content) + "\n")

    print("[=== Contig binning ===] binned into {} clusters".format(im_network.num_top_modules))
    print("[=== Contig binning ===] No. states: {}".format(im_network.num_nodes))
    print("[=== Contig binning ===] No. physical nodes: {}".format(im_network.num_physical_nodes))

    all_nodes_, disputed_ = _pick_the_overlapped()
    fa_content = {}
    for i_header, (i_bin, _) in all_nodes_.items():
        i_seq = contig_fa_lookup[i_header]
        if i_bin in fa_content:
            fa_content[i_bin].append(fasta_concat(i_header, i_seq))
        else:
            fa_content[i_bin] = [fasta_concat(i_header, i_seq)]

    _write_bin_dict(fa_content, out_dir)
    _write_infomap_unbinned_contigs(set(all_nodes_.keys()), contig_fa_lookup, out_dir)
    _write_disputed(disputed_)















