import subprocess
import os

def generate_output(output_prefix, output_format, threads):
    # Execute samtools command to extract hits
    samtools_cmd = f"samtools view -@ {threads} -F 4 -q 5 -h mapped.bam | cut -f1 > mapped.bam_hits.txt"
    try:
        subprocess.run(samtools_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing samtools command: {e}")
        return

    # Rest of the generate_output code...
