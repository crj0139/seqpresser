# Stapler is a tool for mapping SRA reads to reference genomes and quickly appending coordinate info to reads for downstream graphing
# To do - single empty > at beginning of finished reads (line 1)
import subprocess
import os

repodir = os.getcwd()



# INPUT 1 - sample.fastq (or .fasta if no need for conversion)

seqtk_path = "/path/to/seqtk"
input_file = "sample.fastq"
output_file = "sample.fasta"

convert = "seqtk seq -A sample.fastq > sample.fasta" 
subprocess.run(convert, shell=True, cwd=repodir)
