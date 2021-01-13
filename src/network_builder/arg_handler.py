import os
from pathlib import Path
import sys
from shutil import rmtree


def arg_handler(args):
    #     layer_names = ['binned', "assembly", "compo", "pair", "lr"]
    # 1. contig min
    # 2. output path
    # 3.
    # handling contig min length
    ah_min_ctg_len(args)
    # handling output dir

    ah_output(args)
    ah_aln(args)
    ah_read(args)
    args.is_binned, args.is_assembly, args.is_compo = [False] * 3

    ah_bin_dir(args)
    ah_lr_reads(args)
    ah_pairing(args)
    # args.ctg_fa is required
    # handling reads and aln files:

    # args.is_pair and args.is_lr is already set by argparse
    # if args.bin_dir:
    #     args.is_binned = True
    #
    # if args.is_lr:
    #     ah_lr_reads(args)

    return args


def ah_min_ctg_len(args):
    if args.min_ctg_len < 1500:
        args.min_ctg_len = 1500
        print("[::Warning::] minimum contig size must be >= 1500. set it to 1500.")


def ah_output(args):
    if not args.out_dir:
        os.makedirs("netbin_out")
        args.out_dir = os.getcwd() + "/" + "netbin_out/"
    elif os.path.exists(args.out_dir):
        if not args.out_dir.endswith("/"):
            args.out_dir += "/"
        if os.listdir(args.out_dir):
            args.out_dir += "netbin_out/"
            if os.path.exists(args.out_dir):
                rmtree(args.out_dir)
            os.makedirs(args.out_dir)

            print("[::Warning::] directory not empty. will write to or overwrite \"{}\"".format(args.out_dir))
    else:
        if not args.out_dir.endswith("/"):
            args.out_dir += "/"
        Path(args.out_dir).mkdir(parents=True, exist_ok=False)

    Path(args.out_dir + "bins").mkdir(parents=True, exist_ok=False)


def ah_aln(args):
    aln_file = args.aln_file
    if aln_file:
        if len(aln_file) > 2:
            sys.exit("[::Error::] only 1 or 2 sam files are allowed. "
                     "if 2, they should be correspondent to a paired end dataset.")
        elif len(aln_file) == 2:
            args.aln_r1 = aln_file[0]
            args.aln_r2 = aln_file[1]
            del args.aln_file
            print("[::Checking alignment files::] R1:{0}, \n\tR2:{1}".format(args.aln_r1, args.aln_r2))
        elif len(aln_file) == 1:
            args.aln_file = aln_file[0]
            print("[::Checking alignment files::] {}".format(args.aln_file))


def ah_read(args):
    # if paired end in 2 files: args.paired_end_reads does not exist. there is args.reads_r1, args.reads_r2
    # if paired end in 1 file: args.paired_end_reads exist. stub: This is not implemented in the arg parser yet.
    # if single end: there is args.single_read_file
    args.reads_exist = False
    if args.pe_reads:
        reads = args.pe_reads
        args.se_reads = False
        if len(reads) == 1:
            args.reads_exist = True
            # pass
            # stub
            # read paired end reads in interleaved format
        elif len(reads) == 2:
            args.reads_exist = True
            # pass
            # stub
            # read paired end R1 and R2 in separate files
        else:
            print("[::Error::] invalid paired end reads:{}".format(reads))
            del args.pe_reads

    elif args.se_reads:
        args.reads_exist = True
        args.pe_reads = False


def ah_bin_dir(args):
    if args.bin_dir:
        args.is_binned = True
    formatted_bin_dir = []
    for i_bin_dir in args.bin_dir:
        if not i_bin_dir.endswith("/"):
            i_bin_dir += "/"
        formatted_bin_dir.append(i_bin_dir)
    args.bin_dir = formatted_bin_dir


def ah_lr_reads(args):
    # this can only be executed after ah_reads(args)
    if args.is_lr:
        if not args.reads_exist:
            args.is_lr = False
            print("[::Warning::] Read file not provided. Disabling read cloud")
        if not args.aln_file and not args.aln_r1:
            print("[::Calling BWA::]")
        # stub


def ah_pairing(args):
    if args.is_pair:
        if not args.pe_reads and not args.aln_file and not args.aln_r1:
            args.is_pair = False
            print("[Warning] Neither the read file(s) or read to contig alignments provided. Disabling read pairing")
        if args.pe_reads and not args.aln_file and not args.aln_r1:
            print("[Calling BWA]")
            # stub


def bwa_caller():
    # stub
    pass


def samtools_caller():
    # stub
    pass