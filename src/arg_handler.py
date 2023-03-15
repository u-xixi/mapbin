import os
from pathlib import Path
import sys
from shutil import rmtree
from map2bin import *


def arg_handler(args):
    # attributes: fasta file, is_binned, is_assembly, is_pair, min_ctg_len
    #     layer_names = ['binned', "assembly", "compo", "pair", "lr"]
    # 1. contig min
    # 2. output path
    # 3.
    # handling contig min length
    ah_min_ctg_len(args)
    # handling output dir

    ah_output(args)
    ah_out_func(args)
    # ah_read(args)
    args.is_binned, args.is_compo = [False] * 2

    ah_bin_dir(args)
    # ah_lr_reads(args)
    ah_pairing(args)
    ah_assembly(args)
    # args.ctg_fa is required
    # handling reads and aln files:

    # args.is_pair and args.is_lr is already set by argparse
    # if args.bin_dir:
    #     args.is_binned = True
    #
    # if args.is_lr:
    #     ah_lr_reads(args)

    return args

def ah_contigs(args):
    if not args.ctg_fa and not args.bin_dir:
        sys.exit("[=== Error ===] contig files not provided. Exiting...")
#     this function is not used yet
#     because in the arg parser, a contig file is always required



def ah_min_ctg_len(args):
    if args.min_ctg_len < 1000:
        args.min_ctg_len = 1000
        print("[=== Warning ===] minimum contig size must be >= 1000. set it to 1000.")


def ah_output(args):
    if not args.out_dir:
        os.makedirs("mapbin_out")
        args.out_dir = os.getcwd() + "/" + "mapbin_out/"
    elif os.path.exists(args.out_dir):
        if not args.out_dir.endswith("/"):
            args.out_dir += "/"
        if os.listdir(args.out_dir):
            args.out_dir += "mapbin_out/"
            if os.path.exists(args.out_dir):
                rmtree(args.out_dir)
            os.makedirs(args.out_dir)

            print("[=== Warning ===] directory not empty. will write to or overwrite \"{}\"".format(args.out_dir))
    else:
        if not args.out_dir.endswith("/"):
            args.out_dir += "/"
        Path(args.out_dir).mkdir(parents=True, exist_ok=False)

    Path(args.out_dir + "bins").mkdir(parents=True, exist_ok=False)
    Path(args.out_dir + "trivial_bins").mkdir(parents=True, exist_ok=False)
    Path(args.out_dir + "unbinned").mkdir(parents=True, exist_ok=False)
    args.unbinned_short_file = "".join([args.out_dir, "unbinned/unbinned.short.fa"])

    args.init_clu = args.out_dir + "infomap.initial.clu"
    args.im_pajek = args.out_dir + "infomap.net"
    args.im_ftree = args.out_dir + "infomap.ftree"
    args.im_clu = args.out_dir + "infomap.clu"


def ah_out_func(args):
    overlap_type = args.overlap_type
    if overlap_type == 'overlap':
        args.map2bin = map2bin_overlap
    elif overlap_type == 'nodispute':
        args.map2bin = map2bin_nodispute
    elif overlap_type == 'basic':
        args.map2bin = map2bin_basic


def ah_read(args) -> None:
    # not used
    # if paired end in 2 files: args.paired_end_reads does not exist. there is args.reads_r1, args.reads_r2
    # if paired end in 1 file: args.paired_end_reads exist. stub: This is not implemented in the arg parser yet.
    # if single end: there is args.single_read_file
    # args.reads_exist = False
    if args.se_reads:
        print("[=== Info ===] Found single-end reads {}".format(args.se_reads))
    if args.pe_reads:
        print("[=== Info ===] Found paired-end reads {}".format(args.pe_reads))
        if len(args.pe_reads) == 1:
            print("[=== Info ===] Reads as interleaved {}".format(args.pe_reads))
        elif len(args.pe_reads) == 2:
            print("[=== Info ===] Reads as uninterleaved duo {}".format(args.pe_reads))

        else:
            print("[=== Warning ===] discard unexpected paired end reads:{}".format(args.pe_reads))
            args.pe_reads = False


def ah_bin_dir_old(args):
    if args.bin_dir:
        print("[=== Message ===] Found binning result: {}".format(args.bin_dir))
        print("[=== Message ===] \"Binned\" mode on")
        args.is_binned = True
        formatted_bin_dir = []
        for i_bin_dir in args.bin_dir:
            if not i_bin_dir.endswith("/"):
                i_bin_dir += "/"
            formatted_bin_dir.append(i_bin_dir)

        args.bin_dir = formatted_bin_dir


def ah_bin_dir(args):
    if args.bin_dir:
        print("[=== Message ===] Found binning result: {}".format(args.bin_dir))
        print("[=== Message ===] \"Binned\" mode on")
        args.is_binned = True
        if not args.bin_dir.endswith("/"):
            args.bin_dir += "/"


# def ah_lr_reads(args):
#     # this can only be executed after ah_reads(args)
#     if args.is_lr:
#         print("[=== Message ===] Found read cloud file(s): {}".format(args.bin_dir))
#         print("[=== Message ===] \"Read cloud\" mode on")
#         if not args.pe_reads and not args.se_reads:
#             args.is_lr = False
#             print("[=== Warning ===] Read file not provided. Disabling read cloud")
#         if not args.aln_file and not args.aln_r1:
#             print("[=== Calling BWA ===]")
#         # stub

def ah_pairing(args) -> None:
    if args.is_pair:
        if args.aln_files:
            print("[=== Message ===] \"Pairing\" mode on.")
            print("[=== Message ===] Found alignment file(s): {}".format(args.aln_files))


def ah_assembly(args) -> None:
    if args.is_assembly:
        if args.assembler == "spades":
            if args.gfa:
                print("[=== Message ===] \"Assembly\" mode on.")
                print("[=== Message ===] Assembler: Spades. contig gfa file: {}".format(args.gfa))
            else:
                args.is_assembly = False
                print("[=== Message ===] Cannot find contig paths file. Disabling assembly mode...")

