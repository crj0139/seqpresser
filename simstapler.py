# Test on SARS2 https://www.ncbi.nlm.nih.gov/assembly/GCF_009858895.2 SRR14751989
# python simstapler.py -r <reference_fasta>  -t <threads> -o <output_prefix> -d <distribution> -c <coverage>
# python simstapler.py -r SARS2.fna -t 2 -o SARS2testsim -f 1
import argparse
import subprocess

def measure_genome_size(reference_fasta):
    # Run shell command to measure genome size
    command = f"grep -v '>' {reference_fasta} | wc | awk '{{print $3-$1}}'"
    output = subprocess.check_output(command, shell=True, text=True).strip()
    genome_size = int(output)

    return genome_size

def simulate_stapler(reference_fasta, threads, output_prefix, distribution, coverage):
    # Measure genome size
    genome_size = measure_genome_size(reference_fasta)
    print("Genome Size:", genome_size)

    # Run shell command to simulate using seqrequester
    simulate_command = f"seqrequester simulate -genome {reference_fasta} -genomesize {genome_size} -coverage {coverage} -distribution {distribution} > {output_prefix}.sim.fa"
    subprocess.run(simulate_command, shell=True)

    print("Simulation completed!")

def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='SimStapler: A simulation program')

    # Define command-line arguments
    parser.add_argument('-r', '--reference', metavar='<reference_fasta>', required=True, help='Reference FASTA file')
    parser.add_argument('-t', '--threads', metavar='<threads>', type=int, required=True, help='Number of threads')
    parser.add_argument('-o', '--output', metavar='<output_prefix>', required=True, help='Output file prefix')
    parser.add_argument('-d', '--distribution', metavar='<distribution>', required=True, help='Distribution')
    parser.add_argument('-c', '--coverage', metavar='<coverage>', type=int, required=True, help='Coverage')

    # Parse command-line arguments
    args = parser.parse_args()

    # Call the simulation function with the parsed arguments
    simulate_stapler(args.reference, args.threads, args.output, args.distribution, args.coverage)

if __name__ == '__main__':
    main()












#########
#####install seqrequester
#seqrequester simulate -genome <reference_fasta> -genomesize <genome_size> -coverage <coverage> -distribution <read_distribution> > <output_prefix>.sim.fa

#reformat headers