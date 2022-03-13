from arg_handler import *
from network_utils import *
from parsers import *
import argparse
from map2bin import *


def binning(args):
    # step 0: handling args and decide which layer to add
    args = arg_handler(args)

    # step 1: build a network
    # ctg_lookup, unbinned, im_network = generate_ctg_lookup(args.ctg_fa, args.min_ctg_len)
    contig_fa_lookup = parse_contig_fa(args.ctg_fa)
    ctg_lookup, im_network, unbinned_short = generate_ctg_lookup(contig_fa_lookup, args.min_ctg_len)


    # step 2: add links
    if args.is_binned:
        make_infomap_binned(args, im_network=im_network,
                            layer_idx=1, ctg_lookup=ctg_lookup)

    if args.is_assembly:
        make_infomap_assembly(args, im_network,
        layer_idx=2, all_contigs=set(contig_fa_lookup.keys()), ctg_lookup=ctg_lookup)

    if args.is_pair:
        make_infomap_pairing(args, im_network=im_network,
                             layer_idx=3, ctg_lookup=ctg_lookup)

    if args.is_lr:
        make_infomap_linked(args, im_network=im_network,
                            layer_idx=4, ctg_lookup=ctg_lookup)



    # set_im_node_names(ctg_lookup, im_network)
    im_network.run(ftree=True)
    im_network.write_pajek(args.im_pajek)
    im_network.write_flow_tree(args.im_ftree)
    im_network.write_clu(args.im_clu)

    write_unbinned_short(unbinned_short, contig_fa_lookup, args.unbinned_short_file)

    args.net2bin(im_network, contig_fa_lookup, args.out_dir)
    # step 3: run infomap
    # step 4: turn modules into bins

def command_line_parser():
    #     layer_names = ['binned', "assembly", "compo", "pair", "lr"]
    parser = argparse.ArgumentParser()
    #  basics:
    parser.add_argument("-c", "--contigs", required=True, dest="ctg_fa",
                        metavar="contig-fasta", help="assembly contigs in fasta format")
    parser.add_argument("-o", "--out-dir", required=True, dest="out_dir",
                        metavar="out-dir", help="output directory")
    parser.add_argument("-mcl", "--min-ctg-len", type=int, default=2500, dest="min_ctg_len",
                        metavar="contig-min-length", help="minimum contig size. must be >= 1000")
    # output function
    parser.add_argument("-f", "--overlap-type", dest="overlap_type", choices=['overlap', 'nodispute', 'basic'],
                        default='overlap',
                        metavar="contig-overlap", help="whether allowing contig overlap in bins or not.")


    #  integration of binning result. activated by flag -b
    parser.add_argument("-b", "--bins", dest="bin_dir", nargs="+",
                        metavar="binning-dir",
                        help="enable binning integration. Require directory of binning output. can be multiple")
    parser.add_argument("-sfx", "--bin-suffix", default="fa", dest="bin_suffix",
                        metavar="binned-fasta-extension", help="binned fasta file extensions. default: fa")
    parser.add_argument("-csl", "--clique-size-limit", default=25, dest="s", type=int,
                        metavar="clique-size-limit", help="allowed biggest clique for binning network. default: 25")

    #  using assembly graph.
    parser.add_argument("-a", "--assembly", action="store_true", dest="is_assembly",
                        help="integrate assembly graph. Need to specify assembler and relevant graph files.")
    parser.add_argument("-A", "--assembler", dest="assembler", choices=['spades'],
                        help="Assembler used to generate the graph.")
    parser.add_argument("-P", "--contig-paths", dest="ctg_paths",
                        help="spades paths file")
    parser.add_argument("-gfa", "--assembly-graph", dest="gfa",
                        metavar="assembly-gfa",
                        help="enable assembly graph integration. Require assembly graph in gfa format")

    #  using linked reads
    parser.add_argument("-lr", "--linked-reads", action="store_true", dest="is_lr",
                        help="enable integration of read cloud info. Read file required")
    parser.add_argument("-bsl", "--bar-size-limit", type=int, default=2, dest="min_bar_size",
                        help="minimum read cloud size. Default: 2")

    # using read pairing information
    parser.add_argument("-p", "--pairing", action="store_true", dest="is_pair",
                        help="integrate read pairing info. Require paired end reads or their alignments to contigs.")

    #  using read files. for pairs or barcodes.
    #  !!!!fix this: now it accept only 2 sep files. should be able to handle interleaved file too.
    parser.add_argument("-pe", "--paired-end", nargs="+", dest="pe_reads",
                        metavar="paired-reads", help="Still under development. paired end read fastq in 2 separate files")
    parser.add_argument("-se", "--single_end", dest="se_reads",
                        metavar="single-read", help="single-end read fastq")

    #  using sam/bam files. must have if using read clouds or pairing info.
    parser.add_argument("-aln", "--alignment", dest="aln_file",
                        metavar="alignment-sorted-indexed-bam",
                        help="reads to contigs alignments. bam sorted by leftmost coordinates and indexed. one file for both r1 and r2 reads.")


    parser.set_defaults(func=binning)

    return parser.parse_args()


def main():
    args = command_line_parser()
    args.func(args)


if __name__ == '__main__':
    main()