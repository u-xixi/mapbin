from arg_handler import *
from network_utils import *
from parsers import *
import argparse


def binning(args):
    # step 0: handling args and decide which layer to add
    args = arg_handler(args)

    # step 1: build a network
    # ctg_lookup, unbinned, im_network = generate_ctg_lookup(args.ctg_fa, args.min_ctg_len)
    contig_fa_lookup = parse_contig_fa(args.ctg_fa)
    ctg_lookup, unbinned, im_network = generate_ctg_lookup_alt(contig_fa_lookup, args.min_ctg_len)

    # step 2: add links
    if args.is_binned:
        make_infomap_binned(args, im_network=im_network,
                            layer_idx=1, ctg_lookup=ctg_lookup)

    if args.is_pair:
        make_infomap_pairing(args, im_network=im_network,
                             layer_idx=2, ctg_lookup=ctg_lookup)

    if args.is_lr:
        make_infomap_linked(args, im_network=im_network,
                            layer_idx=3, ctg_lookup=ctg_lookup)

    if args.is_assembly:
        make_infomap_assembly(args, im_network,
        layer_idx=4, all_contigs=set(contig_fa_lookup.keys()), ctg_lookup=ctg_lookup)

    # run_infomap(im_network)

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
                        metavar="contig-min-length", help="minimum contig size. must be >= 1500")

    #  integration of binning result. activated by flag -b
    parser.add_argument("-b", "--bins", dest="bin_dir", nargs="+",
                        metavar="binning-dir",
                        help="enable binning integration. Require directory of binning output. can be multiple")
    parser.add_argument("-sfx", "--bin-suffix", default="fa", dest="bin_suffix",
                        metavar="binned-fasta-extension", help="binned fasta file extensions. default: fa")

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
    parser.add_argument("--lr-min-size", type=int, default=2, dest="min_bar_size", help="minimum read cloud size")

    # using read pairing information
    parser.add_argument("-p", "--pairing", action="store_true", dest="is_pair",
                        help="integrate read pairing info. Require paired end reads to contigs alignments. See detail.")

    #  using read files. for pairs or barcodes.
    #  !!!!fix this: now it accept only 2 sep files. should be able to handle interleaved file too.
    parser.add_argument("-pe", "--paired-end", nargs="+", dest="pe_reads",
                        metavar="paired-reads", help="paired end read fastq in 2 separate files")
    parser.add_argument("-se", "--single_end", dest="se_reads",
                        metavar="single-read", help="single-end read fastq")

    #  using sam/bam files. must have if using read clouds or pairing info.
    parser.add_argument("-aln", "--alignment", nargs="+", dest="aln_file",
                        metavar="alignment-sam",
                        help="reads to contigs alignments. sam/bam. one file for all, or 2 files for a pair.")


    parser.set_defaults(func=binning)

    return parser.parse_args()


def main():
    args = command_line_parser()
    args.func(args)


if __name__ == '__main__':
    main()