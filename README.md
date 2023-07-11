# seqpresser
Tools for preparing sequencing data for deep learning

### Dependencies
- [bedtools](https://github.com/arq5x/bedtools2)
- [minimap2](https://github.com/lh3/minimap2)
- [samtools](https://github.com/samtools/samtools)
- [seqrequester](https://github.com/marbl/seqrequester)
- [seqtk](https://github.com/lh3/seqtk)
- [sra-tools](https://github.com/ncbi/sra-tools)



## stapler 
Maps SRA reads to reference genomes, quickly generating alignments for downstream graphing.  Appends, or "staples", alignment information to fasta headers for quick usage.  Also useful for working with public sequencing data in general.

#### Usage

The input for stapler is an SRA accession number, a desired target reference for the SRA reads to be mapped to,
and the desired prefix for output files:

```
python stapler.py -i <input_SRA> -r <reference_fasta>  -t <threads> -o <output_prefix> -f <format_preset>

	-i <input_SRA>: SRA number of reads to be mapped to reference
	-r <reference_fasta>: target reference genome for SRA reads to be mapped to
	-t <threads>: number of threads for mapping
	-o <output_prefix>: Desired prefix for output files
	-f <format_preset>: Preset for desired output format.
		1	>Read_ID strand=<+/->, start=<start_coords>, end=<end_coords> 
```

This program will download the designated SRA data, map it to the target reference, and create a BED file containing
IDs and coordinates.  Furthermore, it will append the mapping coordinate info to the fasta headers, creating a 
".stapled.fasta" output file.


## simstapler
Simulates read sequences from a reference genome using custom distributions. Staples alignment information to simulated reads, and generates a BED file for reference.

#### Usage

```
python simstapler.py -r <reference_fasta>  -t <threads> -o <output_prefix> -d <distribution> -c <coverage>

	-r <reference_fasta>: genome source of simulated reads
	-t <threads>: number of threads for simulating
	-o <output_prefix>: Desired prefix for output files
	-d <distribution>: Values for read length distribution
	-c <coverage>: Simulated coverage of reads
	-f <format_preset>: Preset for determining .fasta header format
		0	>read=<read_id>,<forward/reverse>,position=<start_coords>-<end_coords>,length=<length>,<locus_id>
		1	><read_id> strand=<+/->, start=<start_coords>, end=<end_coords>
```

This will create a BED file with alignment information in it, a set of reads (in .fasta format), and will generate a .fasta file with alignment coordinates stapled to the headers.