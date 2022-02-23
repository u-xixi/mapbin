# Netbin: Metagenomic Contig Binning Using Infomap

## Requirements
- Python 3.8 or 3.9
- [infomap >= 2.0.0](https://github.com/mapequation/infomap) (`pip install infomap`)
- [pysam 0.18.0](https://pysam.readthedocs.io/en/latest/index.html) (`pip install pysam`)
- Lower versions might still work, but no promises

## Usage
The program does 3 main things:
1. build a multilayer network
2. run infomap to partition the network into modules
3. finalize the modules and write out bin

**Example usage:**
```
python netbin/src/network_builder/infomap_binning.py -c contigs.fasta -mcl 3000 \
-b metabat_bins -sfx fa -o netbin_out \
-a -A spades -P contigs.paths \
-p -aln reads2contigs.bam \
-lr -se single_end_reads.fq.gz \
-pe paired_end_r1.fq.gz paired_end_r2.fq.gz
```
You must set `-c <fasta path>` and `-o <out directory>`.

You get to choose what you want to use to build the network of contigs in Step 1. This carries the spirit of a Dominos pizza: the base is the same network clustering algorithm, Infomap, but you could choose your own toppings, use the features of your choice to construct the network. The following features are available, and you need to set them in the command line:
- **Existing binning result** `-b <binning result directory>`
- **An assembly graph** `-a` and you have to set `-A <assembler> -gfa <gfa_file>` or `-P <spades_contig_paths>` for the assembly graph files too
- **read pairing** `-p`, only when your contigs are made with paired-end reads. You need to provide either the reads-to-contigs alignment file, using `-aln <file name>`, or paired-end reads using `-pe <interleaved file>` or `-pe <separated r1 and r2 files>`
- **read cloud linkage** `-lr`, you need to provide reads and alignment files too. Your read file, in fasta or fastq format, should have headers that include a barcode. The barcode should be the last field of the header, delimited by a white space, e.g.: in `ST-J12345:123:ABCDEFGHI:1:1234:12345:67890 ABCDEF`, `ABCDEF` is the barcode. Read more about the linked-read sequencing [here](https://pubmed.ncbi.nlm.nih.gov/26829319/).


**To set the reads:**

You can set paired-end reads using `-pe` or single-end reads using `-se`. You can give `-pe` 1 file to be interpreted as interleaved, or 2 files to be interpreted as R1 and R2 files. `-pe` and `-se` are not exclusive of each other, and can be used at the same time.

Use `python netbin/src/network_builder/infomap_binning.py -h` for more command-line options.


## Output
- `bins/`
which contains the binning result. (Currently only the contigs that are below the length threshold are written out there, as unbinned.short.fa)
- `infomap.ftree` 
Flow tree format, as explained by the Infomap devs,
> The Tree format with an appended section of links with the flow between nodes within the same parent.
- `multilayer.net`
the unpartitioned network in Pajek format. It records the nodes and the links between the nodes. Mainly for debugging. Note that the nodes there are state nodes, see in [here](https://github.com/u-xixi/netbin#get-the-final-binning)
- `infomap.clu` gives you the best assignment of nodes in the top module. It probably can give you the first impression of the network partitioning. This is how it looks:
> \# node_num_id module flow </br>
> 123 7 0.005

## Data for Daniel
Some testing files are available on osa. 
- contig fasta: `/abscratch/xi/haplotagging-May2019/02_assembly/sample1/contigs.fasta`
- Spades contig paths: `/abscratch/xi/haplotagging-May2019/02_assembly/sample1/contigs.paths`
- contig binning results: `/abscratch/xi/haplotagging-May2019/03_binning/sample1/metabat_bins/` or
`/abscratch/xi/haplotagging-May2019/03_binning/sample1/concoct_out/concoct_bins`
- reads to contig alignment: `/abscratch/xi/haplotagging-May2019/04_reads2contigs/sample1.bam`

## Dev status
Currently I have ticked the box for the #1 and #2 functionalities above, and 3 is up for discussion. Read about the details of the implementation in below.
### Building the network
Binning is commonly achieved by analyzing the similarity or associations between contigs, based on which the contigs are clustered into bins. There are a few attributes of contigs that are used for such analysis, including tetranucleotide frequencies, read coverage, contig connectivity (MetaBAT, CONCOCT, MAXBIN) in assembly graphs (GraphBin), read pairing between contigs (Binnacle). There are also a few binning refinement algorithms (Das Tool, GraphBin) which takes the binning result from other programs and try to make some improvement.

In this program, you can choose what features you want to use to cluster contigs. Each feature ends up in a layer. The program build the network layer by layer. This is to get ready for Infomap

Read cloud linkage parts may need some fixes, to work with the newly updated code blocks.
### Running Infomap
The engine is the [Infomap Python API](https://github.com/mapequation/infomap). You should have it installed and the version should better be above 2.0.0.
### Get the final binning
It's a long story...The short version is how it should be done is yet to be decided. The full story is, in a basic network, each objects (nodes) exist in their one and only physical form, and could be assigned to only one module. But in Infomap, the algorithm allows one physical node to have multiple states. Therefore one node, that stands for one contig in our case, can end up in multiple modules. As you may feel the same, people prefers seeing a contig goes to no more than one bin.

A straightforward option is to choose the one with larger flow value in the `.clu` output of infomap. Details in Usage.
