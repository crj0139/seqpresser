import subprocess
import os

def generate_output(output_prefix, output_format):
    # Append coordinates to fasta headers
    mapped_fa_file_path = f"{output_prefix}.mapped.fa"  # Update with the desired filename for the mapped fasta file

    try:
        # Append coordinates to fasta headers using awk command
        append_cmd = (
            f"awk 'NR==FNR{{des[\">\"$1]=$0;next}}/^>/ && des[$1]{{$0=\">\"des[$1]}}1' "
            f"{unique_hit_bed_file_path} {dedup_fasta_file_path} > {mapped_fa_file_path}"
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
            f"sed -e 's/>*-/strand=minus,/g' {mapped_fa_file_path} > {output_prefix}.formatted.fasta"
        )
    elif output_format == 2:
        # Placeholder for format 2 command
        format_cmd = "echo Placeholder for output format 2"
    elif output_format == 3:
        # Placeholder for format 3 command
        format_cmd = "echo Placeholder for output format 3"
    else:
        print("Invalid output format specified.")
        return

    try:
        # Execute the format command
        subprocess.run(format_cmd, shell=True, check=True)
        print("Fasta headers formatted.")
    except subprocess.CalledProcessError as e:
        print(f"Error executing format command: {e}")
        return
