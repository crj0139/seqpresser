import subprocess
import os


def execute_minimap2(reference_genome, fasta_file, threads, output_prefix, output_format):
    # Execute minimap2 command using subprocess
    minimap2_cmd = f"minimap2 -ax map-hifi -t {threads} --sam-hit-only --secondary=no {reference_genome} {fasta_file}"
    sam_file_path = f"{output_prefix}.sam"
    bam_file_path = f"{output_prefix}.bam"

    try:
        subprocess.run(f"{minimap2_cmd} > {sam_file_path}", shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing minimap2 command: {e}")
        return

    # Convert SAM to BAM using samtools
    samtools_cmd = f"samtools view -S -b -@ {threads} {sam_file_path} > {bam_file_path}"
    try:
        subprocess.run(samtools_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing samtools command: {e}")
        return

    # Remove the .sam file
    try:
        os.remove(sam_file_path)
    except OSError as e:
        print(f"Error removing .sam file: {e}")
