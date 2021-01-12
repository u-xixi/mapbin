from Bio import SeqIO
from make_bin_net import make_bin_net
from make_barcode_net import make_barcode_net
import argparse
from arg_handler import *
from file_io import *

# layer_names = ['binned', "assembly", "compo", "pair", "lr"]
# layer_num_lookup

# infomap network format:
#     # *Vertices 5
#     # # node_id name
#     # 1 "i" 10
#     # 2 "j" 10
#     # 3 "k" 10
#     # 4 "l" 1
#     # 5 "m" 1
#     # *Intra
#     # # layer_id node_id node_id weight
#     # 1 1 4 0.8
#     # 1 4 1 1
#     # 1 1 5 0.8
#     # 1 5 1 1
#     # 1 4 5 1
#     # 1 5 4 1
#     # 2 1 2 0.8
#     # 2 2 1 1
#     # 2 1 3 0.8
#     # 2 3 1 1
#     # 2 2 3 1
#     # 2 3 2 1
#


def make_assembly_net(CTG_LOOKUP, layer_id):
    pass


def make_compo_net(args, CTG_LOOKUP, layer_id):
    pass


def make_pair_net(args,CTG_LOOKUP, layer_id):
    pass


def generate_ctg_lookup(ctg_fa, min_ctg_len):
    ctg_lookup = {}
    unbinned = []
    i = 1
    for seq_record in SeqIO.parse(ctg_fa, "fasta"):
        seq_len = len(seq_record)
        if seq_len < min_ctg_len:
           unbinned.append(seq_record)
        else:
            formatted_id = "_".join(seq_record.id.split("_")[:2])
            cur_nbcontig = NbContig(formatted_id, seq_len, str(i))
            ctg_lookup[cur_nbcontig.id] = cur_nbcontig
        i += 1
    return ctg_lookup, unbinned


def layer_numbering(args):
    # this function assigns each layer its correspondent computation method
    # the layer is only referred to by its layer_id
    # the correspondent computational methods are in layer_func
    # the order of the layers is:
    # binning, assembly, tetranucleo + read coverage, read pairing, linked reads

    arg_boo = [args.is_binned, args.is_assembly, args.is_compo, args.is_pair, args.is_lr]
    layer_names = ['binned', "assembly", "compo", "pair", "lr"]
    print(layer_names)
    print(arg_boo)
    layer_func = [make_bin_net, make_assembly_net, make_compo_net, make_pair_net, make_barcode_net]
    layer_num_lookup = []
    layer_ids = deque(range(1, sum(arg_boo) + 1))

    for k, v in zip(arg_boo, layer_func):
        if k:
            i_layer_id = str(layer_ids.popleft())
            layer_num_lookup.append((i_layer_id, v))
    return layer_num_lookup


def write_short_ctg_2unbinned(unbinned, outdir):
    # this is the first load of unbinned fasta: contigs too short
    outfile = outdir + "bins/unbinned.fasta"

    with open(outfile, "w") as f:
        SeqIO.write(unbinned, f, "fasta")


def write_im_vertices(ctg_lookup, net_outfile):
    n_vertices = len(ctg_lookup)
    make_vertex_line = lambda lookup_val: " ".join([lookup_val.num_id, "\"{}\"".format(lookup_val.id)])
    line0 = "*Vertices {}".format(n_vertices)
    vertex_lines = list(map(make_vertex_line, ctg_lookup.values()))
    outstr = line0 + "\n" + "\n".join(vertex_lines) + "\n"
    plain_writer(outstr, net_outfile)


def write_im_layers(args, ctg_lookup, layer_num_lookup, net_outfile):
    outstr = "*Multilayer\n"
    # layer_num_lookup = [(v, k) for k, v in layer_num_lookup.items()]
    # layer_num_lookup.sort()
    print(len(layer_num_lookup))
    for layer_id, func in layer_num_lookup:
        outstr += func(args, ctg_lookup, layer_id)

    with open(net_outfile, 'a') as f_handle:
        f_handle.write(outstr)




# def load_ctgs(ctg_fa, min_ctg_len):
#     # this function reads contigs.fa, generate the CONTIG_LOOKUP and write out the vertices of network file
#     ctg_lookup = {}
#     i = 0
#     for seq_record in SeqIO.parse(ctg_fa, "fasta"):
#         seq_len = len(seq_record)
#         if seq_len >= min_ctg_len:
#             formatted_id = "_".join(seq_record.id.split("_")[:2])






# def assemble_im_network_file(out_file, vertices, layers):
#     out_content = ""

#
#     vertices = list(vertices)
#     vertices.sort()
#     out_content += "*Vertices {}\n".format(str(len(vertices)))
#     out_content += "\n".join(vertices) + "\n"
#
#     out_content += "*Intra\n"
#     out_content += "# layer_id node_id node_id weight\n"
#     layers = list(layers)
#     layers.sort()
#     out_content += "\n".join(layers)
#
#     plain_writer(out_content, out_file)




def command_line_parser():
    #     layer_names = ['binned', "assembly", "compo", "pair", "lr"]
    parser = argparse.ArgumentParser()
    #  basics:
    parser.add_argument("-c", "--contigs", required=True, dest="ctg_fa",
                        metavar="contig-fasta", help="assembly contigs in fasta format")
    parser.add_argument("-o", "--out-dir", required=True, dest="out_dir",
                        metavar="out-dir", help="output directory")
    parser.add_argument("-mcl", "--min-ctg-len", type=int, default=2500, dest="min_ctg_len",
                        metavar="contig-min-length", help="minimum contig size. must be >= 1500")

    #  integration of binning result. activated by flag -b
    parser.add_argument("-b", "--bins", dest="bin_dir", nargs="+",
                        metavar="binning-dir", help="enable binning integration. Require directory of binning output. can be multiple")
    parser.add_argument("-sfx", "--bin-suffix", default="fa", dest="bin_suffix",
                        metavar="binned-fasta-extension", help="binned fasta file extensions. default: fa")

    #  using assembly graph. A stub!!! Activated by gfa
    parser.add_argument("-gfa", "--assembly-graph", dest="gfa",
                        metavar="assembly-gfa", help="enable assembly graph integration. Require assembly graph in gfa format")

    #  using linked reads
    parser.add_argument("-lr", "--linked-reads", action="store_true", dest="is_lr", help="enable integration of read cloud info. Read file required")
    parser.add_argument("--lr-min-size", type=int, default=2, dest="min_bar_size", help="minimum read cloud size")

    # using read pairing information
    parser.add_argument("-p", "--pairing", action="store_true", dest="is_pair", help="integrate read pairing info. Require paired end reads to contigs alignments. See detail.")



    #  using read files. for pairs or barcodes.
    #  !!!!fix this: now it accept only 2 sep files. should be able to handle interleaved file too.
    parser.add_argument("-pe", "--paired-end", nargs=2, dest="pe_reads",
                        metavar="paired-reads", help="paired end read fastq in 2 separate files")
    parser.add_argument("-se", "--single_end", dest="se_reads",
                        metavar="single-read", help="single-end read fastq")




    #  using sam/bam files. must have if using read clouds or pairing info.
    parser.add_argument("-aln", "--alignment", nargs="+", dest="aln_file",
                        metavar="alignment-sam", help="reads to contigs alignments. sam/bam. one file for all, or 2 files for a pair.")




    parser.set_defaults(func=build_layers)



    return parser.parse_args()


def build_layers(args):
    print("++++++ This is NetBin under development +++++++")
    args = arg_handler(args)
    layer_num_lookup = layer_numbering(args)
    network_outfile = args.out_dir + "infomap_input.txt"

    CTG_LOOKUP, unbinned_tooshort = generate_ctg_lookup(args.ctg_fa, args.min_ctg_len)
    write_short_ctg_2unbinned(unbinned_tooshort, args.out_dir)
    write_im_vertices(CTG_LOOKUP, network_outfile)
    write_im_layers(args, CTG_LOOKUP, layer_num_lookup, network_outfile)




# def args_handler(args):
#     formatted_binning_dir = []
#     for i in args.binning_dir:
#         if not i.endswith("/"):
#             i += "/"
#         formatted_binning_dir.append(i)
#     args.binning_dir = formatted_binning_dir
#
#     if args.contig_min < 1500:
#         args.contig_min = 1500
#         print("[Warning] minimum contig size must be >= 1500. set it to 1500.")
#
#     if not args.out_dir:
#         os.makedirs("netbin_out")
#         args.out_dir = os.getcwd() + "/" + "netbin_out/"
#     elif os.path.exists(args.out_dir):
#         if not args.out_dir.endswith("/"):
#             args.out_dir += "/"
#         if os.listdir(args.out_dir):
#             os.makedirs(args.out_dir + "netbin_out")
#             args.out_dir += "netbin_out/"
#             print("[Warning] directory not empty. will write out to \"{}\"".format(args.out_dir))
#     else:
#         if not args.out_dir.endswith("/"):
#             args.out_dir += "/"
#         Path(args.out_dir).mkdir(parents=True, exist_ok=False)
#
#     return args


# def build_layers(args):
#     args = arg_handler(args)
#     CONTIG_LOOKUP = generate_ctg_lookup(args.contig_fa)
#     # print("NODE_29_length_289311_cov_79.10496" in CONTIG_LOOKUP)
#     im_layers = set()
#     im_vertices = set()
#     layer_index = 1
#     for i in args.binning_dir:
#
#         bin_vertices, bin_layer = make_bin_net(i, args.bin_suffix, str(layer_index), CONTIG_LOOKUP)
#         im_vertices.update(bin_vertices)
#         im_layers.update(bin_layer)
#         layer_index += 1
#     network_file = args.out_dir + "multilayer_contig_network.txt"
#     assemble_im_network_file(network_file, im_vertices, im_layers)


# the key --> computation function dictionary:
# layer_names = ['binned', "assembly", "compo", "pair", "lr"]
# key2func_name = {'binned': make_bin_net,
#                  # 'assembly': make_assembly_net,
#                  # 'compo': make_compo_net,
#                  # 'pair': make_pair_net,
#                  'lr': make_barcode_net}


class NbContig():
    # contig is the basic unit of this network
    # use this class as the basis
    # it might help reduce the overhead
    # and allows modifying or adding info to the contig

    def __init__(self, id, length=NotImplemented, num_id=NotImplemented):
        self.id = id
        self.length = length
        self.num_id = num_id


def main():
    args = command_line_parser()
    args.func(args)


if __name__ == '__main__':
    main()
