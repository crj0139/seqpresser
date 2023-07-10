import argparse
import subprocess
import sys


def convert_bam(bam_file_path, threads, output_prefix):
    # Execute samtools and cut commands using subprocess
    cmd = f"samtools view -@ {threads} -F 4 -q 5 -h {bam_file_path} | cut -f1 > {output_prefix}.bam_hits.txt"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing samtools and cut commands: {e}")
        sys.exit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert BAM file")
    parser.add_argument('bam_file', help='BAM file path')
    parser.add_argument('threads', type=int, help='Number of threads')
    parser.add_argument('output_prefix', help='Output prefix')
    args = parser.parse_args()

    convert_bam(args.bam_file, args.threads, args.output_prefix)
