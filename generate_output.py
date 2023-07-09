import subprocess
import os

def generate_output(output_prefix, output_format, threads):
    # Execute samtools command to extract hits from the mapped.bam file
    samtools_cmd = f"samtools view -@ {threads} -F 4 -q 5 -h mapped.bam | cut -f1 > mapped.bam_hits.txt"
    try:
        subprocess.run(samtools_cmd, shell=True, check=True)
        print("Hits extracted from mapped.bam.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing samtools command: {e}")
        return

    # Append coordinates to fasta headers
    mapped_fa_file_path = f"{output_prefix}.mapped.fa"  # Update with the desired filename for the mapped fasta file

    try:
        # Append coordinates to fasta headers using awk command
        append_cmd = (
            f"awk 'NR==FNR{{des[\">\"$1]=$0;next}}/^>/ && des[$1]{{$0=\">\"des[$1]}}1' "
            f"{mapped_bed_file_path} {dedup_fasta_file_path} > {mapped_fa_file_path}"
        )
        subprocess.run(append_cmd, shell=True, check=True)

        print(f"Coordinates appended to fasta headers. Mapped fasta file saved as {mapped_fa_file_path}.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing coordinate appending command: {e}")
        return

    # Format fasta headers based on output format
    if output_format == 1:
        # Format for output format 1
        format_cmd = (
            f"sed -e 's
