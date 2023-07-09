# Test on SARS2 https://www.ncbi.nlm.nih.gov/assembly/GCF_009858895.2 SRR14751989
# stapler.py -i/--input -r/--reference  -t/--threads    -o/--output -f/--format
# python stapler.py -i SRR14751989 -r SARS2.fna -t 2 -o SARS2test -f 1

import argparse
from process import process_sra
from mapper import execute_minimap2

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Stapler program")
    parser.add_argument('-i', '--input', help='Input SRA number')
    parser.add_argument('-f', '--format', help='Output format (1, 2, or 3)')
    parser.add_argument('-o', '--output', help='Output prefix')
    parser.add_argument('-r', '--reference', help='Reference genome file')
    parser.add_argument('-t', '--threads', type=int, help='Number of threads')
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

    # Process the reference genome
    if args.reference:
        process_reference_genome(args.reference)
    else:
        print("Please provide -r/--reference argument.")
        return

    # Execute minimap2
    if args.reference and sra_result and args.threads:
        execute_minimap2(args.reference, sra_result.fasta_file_path, args.threads, args.output, args.format)
    else:
        print("Please provide -r/--reference, -i/--input, -t/--threads arguments.")

if __name__ == '__main__':
    main()


