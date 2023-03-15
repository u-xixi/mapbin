# Mapbin: Metagenomic Contig Binning Using Infomap

## Requirements
- Python 3.8 or 3.9
- [infomap >= 2.0.0](https://github.com/mapequation/infomap) (`pip install infomap`)
- [pysam 0.18.0](https://pysam.readthedocs.io/en/latest/index.html) (`pip install pysam`)
- Lower versions might still work, but no promises

## Usage
The program does 4 main things:
1. Compute the links between input sequences and filter the dubious
2. build a multilayer network based on the links
3. run infomap to partition the network into modules
4. finalize the modules and write out bin

**Example usage:**
```
python infomap_binning.py -c contigs.fasta -mcl 3000 \
-b bins -sfx fa -o mapbin_out \
-a -A spades -P contigs.paths \
-gfa assembly_graph_with_scaffolds.gfa \
-p -aln reads2contigs.bam \
```
```
python infomap_binning.py -c contigs.fasta -mcl 3000 \
-b bins -sfx fa -o mapbin_out \
-p -aln reads2contigs.bam \
```

You must set `-c <fasta path>` and `-o <out directory>`.

You can choose what you want to use to build the network of sequences in Step 1. This carries the spirit of a Dominos pizza: the base is the same network clustering algorithm, Infomap, but you could choose your own toppings, use the features of your choice to construct the network. The following features are available, and you need to set them in the command line:
- **Existing binning result** `-b <binning result directory>`
- **An assembly graph** `-a` and you have to set `-A <assembler> -gfa <gfa_file>` or `-P <spades_contig_paths>` for the assembly graph files too
- **read pairing** `-p`, only when your contigs are made with paired-end reads. You need to provide either the reads-to-contigs alignment file, using `-aln <file name>`.


The two features can be used together, or separately. Use `python infomap_binning.py -h` for more command-line options.


## Output
- `bins/`
contains the main binning result.
- `trivial_bins/`
contains smaller bins that has less than 30,000 bp in total.
- `unbinned/`
contains the unbinned input sequences
- `infomap.ftree` 
Flow tree format, as explained by the Infomap devs,
> The Tree format with an appended section of links with the flow between nodes within the same parent.
- `multilayer.net`
the unpartitioned network in Pajek format. It records the nodes and the links between the nodes. Mainly for debugging. Note that the nodes there are state nodes, see in [here](https://github.com/u-xixi/netbin#get-the-final-binning)
- `infomap.clu` gives you the best assignment of nodes in the top module. It probably can give you the first impression of the network partitioning. This is how it looks:
> \# node_num_id module flow </br>
> 123 7 0.005

