# Test on SARS2 https://www.ncbi.nlm.nih.gov/assembly/GCF_009858895.2 SRR14751989
# python simstapler.py -r <reference_fasta>  -t <threads> -o <output_prefix> -d <distribution> -l <length> -c <coverage> -f <format_preset>
# python simstapler.py -r SARS2.fna -t 2 -o SARS2testsim -f 1 -l 2000-20000 -c 40
import argparse
import subprocess
import os
import re
import shutil

def measure_genome_size(reference_fasta):
    # Run shell command to measure genome size
    command = f"grep -v '>' {reference_fasta} | wc -c"
    output = subprocess.check_output(command, shell=True, text=True).strip()
    genome_size = int(output)

    return genome_size

def simulate_stapler(reference_fasta, threads, output_prefix, distribution, length_range, coverage, format_preset):
    # Measure genome size
    genome_size = measure_genome_size(reference_fasta)
    print("Genome Size:", genome_size)

    # Validate the -d and -l arguments
    if (distribution is None and length_range is None) or (distribution is not None and length_range is not None):
        raise ValueError("Either -d or -l argument should be provided, but not both.")

    # Build the seqrequester command based on the provided arguments
    simulate_command = f"seqrequester simulate -genome {reference_fasta} -genomesize {genome_size} -coverage {coverage}"

    if distribution is not None:
        simulate_command += f" -d {distribution}"
    elif length_range is not None:
        min_length, max_length = parse_length_range(length_range)
        simulate_command += f" -length {min_length}-{max_length}"

    simulate_command += f" > {output_prefix}.sim.fa"

    # Run shell command to simulate using seqrequester
    simulate_process = subprocess.run(simulate_command, shell=True, capture_output=True, text=True)
    if simulate_process.returncode != 0:
        print("Simulation failed.")
        print(simulate_process.stderr)
        return

    if format_preset == '1':
        print("Format preset specified, moving to conversion")
        transform_fasta_file(output_prefix)
    elif format_preset == '0':
        print("No header conversion needed")
        intermediate_files_0 = [f"{output_prefix}.sim.fa"]
        for file in intermediate_files_0:
            if os.path.exists(file):
                new_file = os.path.splitext(file)[0] + ".fasta"
                os.rename(file, new_file)
                
    print("Simulated reads saved")
    
def transform_fasta_file(output_prefix):
    input_file = f"{output_prefix}.sim.fa"
    output_file = f"{output_prefix}.sim.fasta"

    with open(input_file, 'r') as input_handle, open(output_file, 'w') as output_handle:
        for line in input_handle:
            if line.startswith('>'):
                header = line.strip()
                transformed_header = transform_header(header)
                output_handle.write(transformed_header + '\n')
            else:
                output_handle.write(line)

    os.remove(input_file)

def transform_header(header):
    match = re.match(r'>read=(.*),(.+),position=(\d+)-(\d+),length=(\d+),(.+)', header)
    if match:
        read_id = match.group(1)
        forward_reverse = match.group(2)
        start = match.group(3)
        end = match.group(4)
        locus_id = match.group(6)

        strand = '+' if forward_reverse == 'forward' else '-'
        transformed_header = f'>{read_id} strand={strand}, start={start}, end={end}'
        return transformed_header

    return header

def parse_length_range(length_range):
    min_length, max_length = map(int, length_range.split('-'))
    return min_length, max_length

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='SimStapler: A simulation program')

    # Define command-line arguments
    parser.add_argument('-r', '--reference', metavar='<reference_fasta>', required=True, help='Reference FASTA file')
    parser.add_argument('-t', '--threads', metavar='<threads>', type=int, required=True, help='Number of threads')
    parser.add_argument('-o', '--output', metavar='<output_prefix>', required=True, help='Output file prefix')
    parser.add_argument('-d', '--distribution', metavar='<distribution>', help='Distribution')
    parser.add_argument('-l', '--length', metavar='<length_range>', help='Length range (min-max)')
    parser.add_argument('-c', '--coverage', metavar='<coverage>', type=int, required=True, help='Coverage')
    parser.add_argument('-f', '--format', metavar='<format_preset>', choices=['0', '1'], default='0', help='Format preset (0 or 1)')

    # Parse command-line arguments
    args = parser.parse_args()

    # Split the length range into min and max values
    length_range = None
    if args.length is not None:
        length_range = args.length

    # Call the simulation function with the parsed arguments
    simulate_stapler(args.reference, args.threads, args.output, args.distribution, length_range, args.coverage, args.format)

if __name__ == '__main__':
    main()
