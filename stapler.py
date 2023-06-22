# Stapler is a tool for mapping SRA reads to reference genomes and quickly appending coordinate info to reads for downstream graphing
# To do - single empty > at beginning of finished reads (line 1)
import subprocess
import os

repodir = os.getcwd()

seqtk = "seqtk seq -A sample.fastq > sample.fasta" 
subprocess.run(command, shell=True, cwd=repodir)
