# Test on SARS2 https://www.ncbi.nlm.nih.gov/assembly/GCF_009858895.2 SRR14751989
# stapler.py -i/--input -r/--reference  -t/--threads    -o/--output -f/--format
# python stapler.py -i SRR14751989 -r SARS2.fna -t 2 -o SARS2test -f 1

import argparse
import subprocess
import os

class SRAImportResult:
    def __init__(self, fastq_file_path, fasta_file_path):
        self.fastq_file_path = fastq_file_path
        self.fasta_file_path = fasta_file_path

def process_sra(sra_number):
    # Execute prefetch command using subprocess
    prefetch_cmd = f"prefetch {sra_number}"
    try:
        subprocess.run(prefetch_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing prefetch command: {e}")
        return

    # Extract the directory name
    directory_name = sra_number

    # Check if the directory exists
    if os.path.exists(directory_name):
        # Change working directory to the prefetch-created directory
        os.chdir(directory_name)

        # Execute fasterq-dump command using subprocess
        fasterq_cmd = f"fasterq-dump {sra_number}"
        try:
            subprocess.run(fasterq_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing fasterq-dump command: {e}")
            return

        # Capture the generated .fastq file path
        fastq_file = f"{sra_number}.fastq"
        fastq_file_path = os.path.join(os.getcwd(), fastq_file)

        # Execute seqtk command to convert .fastq to .fasta
        fasta_file = f"{sra_number}.fasta"
        fasta_file_path = os.path.join(os.getcwd(), fasta_file)
        seqtk_cmd = f"seqtk seq -A {fastq_file_path} > {fasta_file_path}"
        try:
            subprocess.run(seqtk_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing seqtk command: {e}")
            return

        # Change back to the original working directory
        os.chdir('..')

        # Return the SRA import result
        return SRAImportResult(fastq_file_path, fasta_file_path)
    else:
        print(f"Directory '{directory_name}' does not exist.")
        return


def execute_minimap2(reference_genome, fasta_file, threads, output_prefix, output_format):
    # Execute minimap2 command using subprocess
    minimap2_cmd = f"minimap2 -ax map-hifi -t {threads} --sam-hit-only --secondary=no {reference_genome} {fasta_file}"
    sam_file_path = "mapped.sam"
    bam_file_path = "mapped.bam"

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

    # Generate output
    generate_output(output_prefix, output_format, threads)


def generate_output(output_prefix, output_format, threads):
    # Execute samtools command to extract hits
    samtools_cmd = f"samtools view -@ {threads} -F 4 -q 5 -h mapped.bam | cut -f1 > mapped.bam_hits.txt"
    try:
        subprocess.run(samtools_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing samtools command: {e}")
        return

    # Rest of the generate_output code...


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Stapler program")
    parser.add_argument('-i', '--input', help='Input SRA number')
    parser.add_argument('-f', '--format', help='Output format (1, 2, or 3)')
    parser.add_argument('-o', '--output', help='Output prefix')
    parser.add_argument('-r', '--reference', required=True, help='Reference genome file')
    parser.add_argument('-t', '--threads', type=int, default=2, help='Number of threads')
    args = parser.parse_args()

    # Process the input based on the provided arguments
    if args.input:
        sra_result = process_sra(args.input)
        if not sra_result:
            print("Error processing SRA.")
            return
    else:
        print("Please provide -i/--input argument.")
        return

    # Execute minimap2
    if sra_result:
        execute_minimap2(args.reference, sra_result.fasta_file_path, args.threads, args.output, args.format)
    else:
        print("Error processing SRA.")


if __name__ == '__main__':
    main()

