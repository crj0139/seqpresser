<p align="left">
  <img src="misc/seqp.jpg" width="500" title="seqp_logo">
</p>
Seqpresser is a toolkit for preparing sequencing data for deep learning, specifically for downstream neural networks.  
The stapler.py program was designed to be a useful, simple tool for quickly increasing public sequencing read support 
for training graph neural networks that are intended to (re)assemble user-specified reads.  The simstapler.py program is 
intended to quickly generate large amounts of simulated reads from a reference genome(s), which can later be independently 
re-assembled via an OLC-based assembler and compared to the original simulated read source for graph training.

### Dependencies
- [sra-tools](https://github.com/ncbi/sra-tools) (needs SRA  prefetch and fasterq-dump in PATH)
- [bedtools](https://github.com/arq5x/bedtools2)
- [minimap2](https://github.com/lh3/minimap2) (>= 2.26)
- [samtools](https://github.com/samtools/samtools) (>= 1.15.1)
- [seqrequester](https://github.com/marbl/seqrequester)
- [seqtk](https://github.com/lh3/seqtk)

All dependencies except SRA-tools can be installed and loaded via:
```
bash ./install.sh
conda activate seqpresser
```

## stapler 
Maps long SRA reads to reference genomes, quickly generating alignments for downstream graphing.  Appends, or "staples", alignment information to fasta headers for quick reading.
Useful for high-throughput mapping of public sequencing data to a user-specified reference assembly.  Output is designed to
be used for downstream training/validation in assembly graphs.

(Also useful for working with public sequencing data in general).

#### Usage

The input for stapler is an SRA accession number, a desired target reference for the SRA reads to be mapped to,
and the desired prefix for output files:

```
python stapler.py -i <input_SRA> -r ./<reference_fasta>  -t <threads> -o <output_prefix> -f <format_preset>

	-i <input_SRA>: SRA number of reads to be mapped to reference
	-r <reference_fasta>: target reference genome for SRA reads to be mapped to
	-t <threads>: number of threads for mapping
	-o <output_prefix>: Desired prefix for output files
	-f <format_preset>: Preset for desired output format.
		1	>Read_ID strand=<+/->, start=<start_coords>, end=<end_coords> 
	-p <minimap2_preset>: preset for mapping parameters (designated by -ax in minimap2)
```

This program will download the designated SRA data, map it to the target reference, and create a BED file containing
IDs and coordinates.  Furthermore, it will append the mapping coordinate info to the fasta headers, creating a 
".stapled.fasta" output file.


## simstapler
Simulates read sequences from a reference genome using custom distributions. Staples alignment information to simulated reads, and generates a BED file for reference.

#### Usage

```
python python simstapler.py -r ./<reference_fasta>  -t <threads> -o <output_prefix> -d <distribution> -l <length> -c <coverage> -f <format_preset>

	-r <reference_fasta>: genome source of simulated reads
	-t <threads>: number of threads for simulating
	-o <output_prefix>: Desired prefix for output files
	-d <distribution>: Values for read length distribution (only if no -l given)
	-l <length>: Range of read lengths for simulated reads <min-max> (only if no -d given)
	-c <coverage>: Simulated coverage of reads
	-f <format_preset>: Preset for determining .fasta header format
		0	>read=<read_id>,<forward/reverse>,position=<start_coords>-<end_coords>,length=<length>,<locus_id>
		1	><read_id> strand=<+/->, start=<start_coords>, end=<end_coords>
```

This will create a set of reads (in .fasta format) with alignment coordinates stapled to the headers.



