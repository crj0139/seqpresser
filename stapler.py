import os
import argparse
import subprocess

# stapler.py -i/--input -r/--reference  -t/--threads    -o/--output -f/--format

path_variable = os.environ.get('PATH')
print(f"Current PATH variable: {path_variable}")


def process_sra(sra_number):
    # Execute prefetch command using subprocess
    cmd = f"prefetch {sra_number}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing prefetch command: {e}")
        return

    # Continue with the remaining processing logic
    pass








# ################################################

def process_input_genome(input_genome):
    # Placeholder for input genome processing logic
    print(f"Processing input genome: {input_genome}")

def generate_output(output_prefix, output_format):
    # Placeholder for output file generation logic
    print(f"Generating output: {output_prefix}.{output_format}")

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Stapler program")
    parser.add_argument('-s', '--sra', help='SRA number')
    parser.add_argument('-i', '--input', help='Input genome file')
    parser.add_argument('-f', '--format', help='Output format (1, 2, or 3)')
    parser.add_argument('-o', '--output', help='Output prefix')
    args = parser.parse_args()

    # Process the input based on the provided arguments
    if args.sra:
        process_sra(args.sra)
    elif args.input:
        process_input_genome(args.input)
    else:
        print("Please provide either -s/--sra or -i/--input argument.")

    # Generate the output file
    if args.output and args.format:
        generate_output(args.output, args.format)
    else:
        print("Please provide both -f/--format and -o/--output arguments.")

if __name__ == '__main__':
    main()
