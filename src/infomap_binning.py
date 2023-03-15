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
    ctg_lookup, unbinned_short = generate_ctg_lookup(contig_fa_lookup, args.min_ctg_len)
    max_ctg_len = max([v.length for v in ctg_lookup.values()])
    # ctg_numid_lookup = {v.num_id:v for k, v in ctg_lookup.items()}
    layer = {'bin': 1, 'pair': 2, 'assembly': 3, 'linked': 4}

    # # testing pairing
    # if args.is_pair_test:
    #     make_pairing_flat_net(args, ctg_lookup)
    # else:
    im_multi = None
    # step 2: add links
    # if args.bin_as_init:
    #     if args.is_binned:
    #         make_infomap_initial(args, ctg_lookup)
    #     else:
    #         print("[=== Warning ===] No bins provided. ignore --bin-as-initial")

    im_multi = init_infomap(args, ctg_lookup)

    bin_lookup = {}
    if args.is_binned:
        bin_lookup = make_binned_lookup(args, ctg_lookup)
    else:
        bin_lookup = {ctg: 1 for ctg in ctg_lookup}

    if args.is_pair:
        make_infomap_pair(args, im_multi, layer, ctg_lookup, bin_lookup, args.nproc, max_ctg_len)

    if args.is_assembly:
        make_infomap_assembly(args, im_multi, layer,
                              ctg_lookup, bin_lookup, max_ctg_len)

    if args.is_binned:
        make_infomap_binned(im_multi, layer, ctg_lookup, bin_lookup)


    im_multi.run(ftree=True)
    # im_multi.write_pajek(args.im_pajek)
    im_multi.write_state_network(args.im_pajek)
    im_multi.write_flow_tree(args.im_ftree)
    im_multi.write_clu(args.im_clu, states=True)

    # step 4: turn modules into bins
    write_unbinned_short(unbinned_short, contig_fa_lookup, args.unbinned_short_file)
    args.map2bin(im_multi, contig_fa_lookup, args.min_single_len, args.out_dir, args.min_bin_bp)


def command_line_parser():
    #     layer_names = ['binned', "assembly", "compo", "pair", "lr"]
    parser = argparse.ArgumentParser()
    #  basics:
    parser.add_argument("-c", "--contigs", required=True,
                        dest="ctg_fa",
                        metavar="contig-fasta",
                        help="Assembly contigs in fasta format")
    parser.add_argument("-o", "--out-dir", required=True, dest="out_dir",
                        metavar="out-dir", help="Output directory")
    parser.add_argument("-mcl", "--min-ctg-len",
                        type=int, default=1500,
                        dest="min_ctg_len",
                        metavar="contig-min-length",
                        help="Minimum contig size.")
    parser.add_argument("-msl", "--min-single-len",
                        type=int, default=100000,
                        dest="min_single_len",
                        metavar="min-sbin-length",
                        help="Contig above this length can form a single bin. default: 100000.")
    parser.add_argument("-mb", "--min-bin-bp", type=int,
                        dest="min_bin_bp", default=30000,
                        help="bins that has less bp than this will not be in the output.")
    # output function

    parser.add_argument("-f", "--overlap-type",
                        dest="overlap_type",
                        choices=['overlap', 'nodispute', 'basic'],
                        default='overlap',
                        metavar="contig-overlap",
                        help="Whether allowing contig overlap in bins or not.")
    parser.add_argument("-si", "--seed-infomap",
                        default=1, type=int,
                        dest="im_seed",
                        metavar="random-seed-infomap",
                        help="Set seed for running infomap. default: 1")


    #  integration of binning result. activated by flag -b
    # parser.add_argument("-b", "--bins", dest="bin_dir", nargs="+",
    #                     metavar="binning-dir",
    #                     help="enable binning integration. Require directory of binning output. can be multiple")
    parser.add_argument("-b", "--bins", dest="bin_dir",
                        metavar="binning-dir",
                        help="Use precomputed binning result. Require directory of binned fasta.")
    # parser.add_argument("-B", "--bin-as-initial", action="store_true",
    #                     dest="bin_as_init",
    #                     help="Flag to set the binning result as initial partitioning to further optimize.")
    parser.add_argument("-sfx", "--bin-suffix", default="fa",
                        dest="bin_suffix",
                        metavar="binned-fasta-extension",
                        help="binned fasta file extensions. default: fa")

    #  using assembly graph.
    parser.add_argument("-a", "--assembly", action="store_true",
                        dest="is_assembly",
                        help="Integrate assembly graph. "
                             "Need to specify assembler and relevant graph files.")
    parser.add_argument("-A", "--assembler", dest="assembler",
                        choices=['spades'],
                        help="Assembler used to generate the graph.")
    parser.add_argument("-P", "--contig-paths", dest="ctg_paths",
                        help="Spades paths file")
    parser.add_argument("-gfa", "--assembly-graph", dest="gfa",
                        metavar="assembly-gfa",
                        help="Enable assembly graph integration. "
                             "Require assembly graph in gfa format")
    parser.add_argument("-l", "--min-percent-link", type=float,
                        dest="l_perc", metavar="min-percent-link",
                        help="two contigs can only be connected when "
                             "the connection part is larger than this "
                             "percentage of the shorter contig",
                        default=0.1)

    #  using linked reads
    # parser.add_argument("-lr", "--linked-reads", action="store_true", dest="is_lr",
    #                     help="enable integration of read cloud info. Read file required")
    # parser.add_argument("-bsl", "--bar-size-limit", type=int, default=2, dest="min_bar_size",
    #                     help="minimum read cloud size. Default: 2")

    # using read pairing information
    parser.add_argument("-p", "--pairing", action="store_true", dest="is_pair",
                        help="integrate read pairing info. Require paired end reads or their alignments to contigs.")
    # parser.add_argument("-pt", "--pairing-test",
    #                     action="store_true",dest="is_pair_test",
    #                     help="Test pairing network partitioned by 2level algo")
    parser.add_argument("-t", "--ncpus", type=int, dest="nproc",
                        help="Number of CPUs to use.", default=1)

    #  using sam/bam files. must have if using read clouds or pairing info.
    parser.add_argument("-aln", "--alignment", dest="aln_files",
                        nargs="*",
                        metavar="reads-to-contigs",
                        help="Reads to contigs alignments. bam "
                             "sorted by leftmost coordinates and indexed. "
                             "one file for both r1 and r2 reads.")
    parser.add_argument("-aqp", "--aligned-query-percent", dest="aln_query_perc",
                        metavar="aligned-query-percent", default=0.95, type=float,
                        help="aligned-bases/total-bases in query. Default: 0.95")

    parser.add_argument("-mis", "--max-insert-size", dest="max_insert_size",
                        metavar="max-insert-size", default=800, type=int,
                        help="Maximum insert size allowed. Default: 800")

    parser.add_argument("-rl", "--rlen", dest="rlen",
                        metavar="read-length", default=150, type=int,
                        help="Read length. Default: 150")

    parser.set_defaults(func=binning)

    return parser.parse_args()


def main():
    args = command_line_parser()
    args.func(args)


if __name__ == '__main__':
    main()