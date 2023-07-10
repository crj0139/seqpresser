# Test on SARS2 https://www.ncbi.nlm.nih.gov/assembly/GCF_009858895.2 SRR14751989
# stapler.py -i/--input -r/--reference  -t/--threads    -o/--output -f/--format
# stapler.py -i <input_SRA> -r <reference_fasta>  -t <threads>    -o <output_prefix> -f <format_preset>
# python stapler.py -i SRR14751989 -r SARS2.fna -t 2 -o SARS2test -f 1

import argparse
import subprocess
import sys
from process import process_sra
from map import execute_minimap2


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

        # Execute minimap2
        execute_minimap2(
            args.reference,
            sra_result.fasta_file_path,
            args.threads,
            args.output,
            args.format
        )

        # Pass the BAM file path to the conversion script
        if args.format:
            bam_file_path = f"{args.output}.bam"
            convert_cmd = f"python convert.py {bam_file_path} {args.threads} {args.output}"
            try:
                subprocess.run(convert_cmd, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error executing convert.py script: {e}")
                sys.exit(1)

    else:
        print("Please provide -i/--input argument.")

    print("Stapler program completed successfully.")


if __name__ == '__main__':
    main()

