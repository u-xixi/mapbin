import glob
from file_io import *
import sys
from itertools import combinations
import time

#  generate infomap layer using binning results
## input: a dir of contig bins. output: text of infomap layer
## im is short for infomap


# def make_bin_net(bin_dirs, suffix, min_contig_len, layer_id, CTG_LOOKUP):
def make_bin_net(args, CTG_LOOKUP, layer_id):
    bin_dir = args.bin_dir
    suffix = args.bin_suffix
    min_ctg_len = args.min_ctg_len

    def update_dict_val(base_dict, update_dict):
        for key, val in update_dict.items():
            if key in base_dict:
                base_dict[key] += val
            else:
                base_dict[key] = val
    print("[Binning layer] Computing network from binning result.")
    im_edges_content = {}

    for i_bin_dir in bin_dir:
        bin_fa_files = _get_fasta_files_in_dir(i_bin_dir, suffix)
        if not len(bin_fa_files):
            sys.exit("[::Error::] No contigs detected in \"{}\"".format(i_bin_dir))
        i_edges = {}

        time_begin = time.time()
        for bin_fa in bin_fa_files:
            ii_edges = _get_bin_info_content(bin_fa, min_ctg_len, CTG_LOOKUP)
            if set(ii_edges).intersection(i_edges):
                sys.exit("Binning result invalid: bins overlapping.\nRepeated contig detected in file {}".format(bin_fa))
            else:
                i_edges.update(ii_edges)
        if not i_edges:
            print("Warning: no valid contigs detected. Please check the binning result.")
        update_dict_val(im_edges_content, i_edges)
        time_finish = time.time()
        print("[Binning layer] parsing binning result takes {:.1f} secs".format(time_finish - time_begin))

    if im_edges_content:
        out_content = [" ".join([layer_id, k, "{:.5f}".format(v)]) for k, v in im_edges_content.items()]
        out_content.sort()
        return "\n".join(out_content) + "\n"
    else:
        print("Warning: no valid binning result detected. Discard binning result layer.")
        return ""


# def _get_bin_info_content(bin_fa, min_contig_len, layer_id, CONTIG_LOOKUP):
#     # returns a set
#     fa_content = _get_seq_file_content(bin_fa)
#     fa_headers = [_parse_seq_name(i) for i in fa_content if i.startswith(">")]  # headers don't have ">"
#
#     # im_vertices = set([_make_vertex_line(i, CONTIG_LOOKUP) for i in fa_headers])
#     im_edges = set()
#
#     for a, b in combinations(fa_headers, 2):
#         contig_a = CONTIG_LOOKUP[a]
#         contig_b = CONTIG_LOOKUP[b]
#         a2b_line = " ".join([layer_id, contig_a.num_id, contig_b.num_id, "{:5f}".format(min_contig_len/contig_a.length)])
#         b2a_line = " ".join([layer_id, contig_b.num_id, contig_a.num_id, "{:5f}".format(min_contig_len/contig_b.length)])
#         im_edges.add(a2b_line)
#         im_edges.add(b2a_line)
#
#     return im_edges
def _get_bin_info_content(bin_fa, min_contig_len, CTG_LOOKUP):
    # returns a set
    fa_content = _get_seq_file_content(bin_fa)
    fa_headers = [_parse_seq_name(i) for i in fa_content if i.startswith(">")]  # headers don't have ">"
    fa_headers = [i for i in fa_headers if i in CTG_LOOKUP]
    im_edges = {}

    for a, b in combinations(fa_headers, 2):
        ctg_a = CTG_LOOKUP[a]
        ctg_b = CTG_LOOKUP[b]
        a2b_key = " ".join([ctg_a.num_id, ctg_b.num_id])
        a2b_val = min_contig_len/ctg_a.length
        b2a_key = " ".join([ctg_b.num_id, ctg_a.num_id])
        b2a_val = min_contig_len/ctg_b.length
        im_edges[a2b_key] = a2b_val
        im_edges[b2a_key] = b2a_val

    return im_edges


def _get_seq_file_content(file_path):

    seq_file = SeqFile(file_path)
    if seq_file.seqformat == "fasta":
        return seq_file.file_content
    else:
        sys.exit("Invalid Fasta file: {}".format(file_path))


def _get_fasta_files_in_dir(dir, suffix):
    # dir must ends up with "/"
    return glob.glob(dir + "*."+ suffix)


_parse_seq_name = lambda x: "_".join(x[1:-1].split("_")[:2])




