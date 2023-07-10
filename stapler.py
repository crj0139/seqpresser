# Test on SARS2 https://www.ncbi.nlm.nih.gov/assembly/GCF_009858895.2 SRR14751989
# stapler.py -i/--input -r/--reference  -t/--threads    -o/--output -f/--format
# python stapler.py -i <input_SRA> -r <reference_fasta>  -t <threads>    -o <output_prefix> -f <format_preset>
# python stapler.py -i SRR14751989 -r SARS2.fna -t 2 -o SARS2test -f 1
import argparse
import subprocess
import sys
import os
from process import process_sra
from map import execute_minimap2
from convert import convert_bam

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

        # Check if the SRA files already exist
        if sra_result.fasta_file_path:
            print("SRA files already present, skipping processing")

        # Check if the output BAM file already exists
        bam_file_path = f"{args.output}.bam"
        if os.path.exists(bam_file_path):
            print("Alignment .bam already found, skipping mapping")
        else:
            # Execute minimap2
            execute_minimap2(
                args.reference,
                sra_result.fasta_file_path,
                args.threads,
                args.output,
                args.format
            )

        # Check if the output BAM file exists after mapping
        if os.path.exists(bam_file_path):
            # Pass the BAM file path and input SRA to the conversion script
            if os.path.exists(f"{args.output}.uniqhit_dedup.fa"):
                print("Conversion skipped as .bam file already exists.")
            else:
                convert_bam(
                    bam_file_path,
                    args.threads,
                    args.output,
                    args.reference,
                    args.input
                )
        else:
            print("Mapping step failed.")
            return
            
        # Check if the unique reads have been extracted
#        if os.path.exists(reads_uniq_dedup_fa_file_path):
#            # Pass the BAM file path and input SRA to the conversion script
#            if os.path.exists(f"{args.output}.uniqhit_dedup.fa"):
#                print("C")
#            else:

#        else:
#            print("Conversion step failed.")
#            return

    else:
        print("Please provide -i/--input argument.")

    print("Stapler program completed successfully.")


if __name__ == '__main__':
    main()

