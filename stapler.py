
# python stapler.py -i <input_SRA> -r <reference_fasta>  -t <threads>    -o <output_prefix> -f <format_preset>  -p <minimap2_preset>
# python stapler.py -i SRR16925101	 -r ./seqpr_test.fna -t 2 -o seqprtest -f 1 -p map-pb
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
    parser.add_argument('-p', '--minimap2_preset', help='Preset value for minimap2')
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
            # Execute minimap2 with the preset value
            execute_minimap2(
                args.reference,
                sra_result.fasta_file_path,
                args.threads,
                args.output,
                args.format,
                args.minimap2_preset
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

