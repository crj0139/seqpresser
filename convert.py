import argparse
import subprocess
import sys
import os
import glob


def convert_bam(bam_file_path, threads, output_prefix, reference_fasta, input_sra):
    # Execute samtools and cut commands using subprocess
    cmd = f"samtools view -@ {threads} -F 4 -q 5 -h {bam_file_path} | cut -f1 > {output_prefix}.bam_hits.txt"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing samtools and cut commands: {e}")
        return

    # Convert BAM to BED
    bed_file_path = f"{output_prefix}.hits.bed"
    cmd = f"bamToBed -i {bam_file_path} > {bed_file_path}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing bamToBed command: {e}")
        return

    # Format the BED headers
    hits_a_bed_file_path = f"{output_prefix}.hits.a.bed"
    hits_b_bed_file_path = f"{output_prefix}.hits.b.bed"
    hits_c_bed_file_path = f"{output_prefix}.hits.c.bed"
    cmd = (
        f"cut -f 2- {bed_file_path} > {hits_a_bed_file_path} && "
        f"awk 'BEGIN {{FS=OFS=\"\\t\"}} {{print $3, $1, $2, $4, $5}}' {hits_a_bed_file_path} > {hits_b_bed_file_path} && "
        f"awk 'BEGIN {{FS=OFS=\"\\t\"}} {{print $1, $5, $2, $3}}' {hits_b_bed_file_path} > {hits_c_bed_file_path}"
    )
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error formatting BED headers: {e}")
        return

    # Extract uniquely mapped read IDs from BED
    uniqhit_bed_file_path = f"{output_prefix}.uniqhit.bed"
    cmd = f"grep -wFf {output_prefix}.bam_hits.txt {output_prefix}.hits.c.bed > {uniqhit_bed_file_path}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error extracting uniquely mapped read IDs: {e}")
        return

    # Extract sequences using seqtk
    reads_uniq_fa_file_path = f"{output_prefix}.reads_uniq.fa"
    cmd = f"seqtk subseq {input_sra}/{input_sra}.fasta {uniqhit_bed_file_path} > {reads_uniq_fa_file_path}"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error extracting sequences using seqtk: {e}")
        return

    # Remove duplicates from the reads_uniq.fa file
    reads_uniq_fa_file_path = f"{output_prefix}.reads_uniq.fa"
    reads_uniq_dedup_fa_file_path = f"{output_prefix}.reads_uniq_dedup.fa"

    cmd = f"""
    sed -e '/^>/s/$/@/' -e 's/^>/#/' {reads_uniq_fa_file_path} | \
    tr -d '\n' | tr "#" "\\n" | tr "@" "\\t" | \
    sort -u -t ' ' -f -k1,1 | \
    sed -e 's/^/>/' -e 's/\\t/\\n/' > {reads_uniq_dedup_fa_file_path}
    """

    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error removing duplicates from reads_uniq.fa file: {e}")
        return

    # Fix fasta headers and annotate with strand, start, and end positions
    mapped_a_fa_file_path = f"{output_prefix}.mapped.a.fa"
    cmd = (
        f"awk 'NR==FNR{{des[\">\"$1]=$0;next}}/^>/&&des[$1]{{$0=\">\"des[$1]}}1' "
        f"{uniqhit_bed_file_path} {reads_uniq_dedup_fa_file_path} > {mapped_a_fa_file_path}"
    )
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error fixing fasta headers: {e}")
        return

    # Format the header information
    mapped_b_fa_file_path = f"{output_prefix}.mapped.b.fa"
    mapped_c_fa_file_path = f"{output_prefix}.mapped.c.fa"
    mapped_d_fa_file_path = f"{output_prefix}.mapped.d.fa"
    mapped_e_fa_file_path = f"{output_prefix}.mapped.e.fa"
    mapped_f_fa_file_path = f"{output_prefix}.mapped.f.fa"
    mapped_g_fa_file_path = f"{output_prefix}.mapped.g.fa"
    cmd = (
        f"sed -e 's/>*-/strand=minus,/g' {mapped_a_fa_file_path} > {mapped_b_fa_file_path} && "
        f"sed -e 's/>*+/strand=+,/g' {mapped_b_fa_file_path} > {mapped_c_fa_file_path} && "
        f"sed -e 's/>*strand=+,\t/strand=+,\tstart=/g' {mapped_c_fa_file_path} > {mapped_d_fa_file_path} && "
        f"sed -e 's/>*strand=minus,\t/strand=minus,\tstart=/g' {mapped_d_fa_file_path} > {mapped_e_fa_file_path} && "
        f"sed 's/\\(.*\\)\\t/\\1,\\tend=/' {mapped_e_fa_file_path} > {mapped_f_fa_file_path} && "
        f"sed -e 's/>*=minus/=-/g' {mapped_f_fa_file_path} > {mapped_g_fa_file_path}"
    )
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error formatting header information: {e}")
        return
    
    # Fix empty >
    cmd = f"awk 'BEGIN {{RS = \">\" ; FS = \"\\n\" ; ORS = \"\"}} $2 {{print \">\"$0}}' {mapped_g_fa_file_path} > {output_prefix}.int.fa && \
        mv {output_prefix}.int.fa {output_prefix}.blanked.fasta"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error fixing empty > lines: {e}")
        return

    # Remove consecutive whitespaces and blank lines
    cmd = f"sed -e \"s/>*[[:space:]]\+/ /g\" ./{output_prefix}.blanked.fasta > ./{output_prefix}.0.fa && \
        sed -e '/^[[:space:]]*$/d' ./{output_prefix}.0.fa > ./{output_prefix}.stapled.fasta"
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error removing consecutive whitespaces and blank lines: {e}")
        return

    # List of intermediate files to remove
    file_extensions = ["reads_uniq_dedup.fa", "csv", "a.fa", "b.fa", "c.fa", "d.fa", "e.fa", "f.fa", "g.fa",
    "a.bed", "b.bed", "c.bed", "hits.bed", "blanked.fasta", "reads_uniq.fa", "bam_hits.txt", "0.fa"]

    # Remove files with specified extensions (excluding output_prefix files)
    for file_extension in file_extensions:
        files_to_remove = glob.glob(f"./*.{file_extension}")
        for file_path in files_to_remove:
            try:
                os.remove(file_path)
                print(f"Removed file: {file_path}")
            except FileNotFoundError:
                pass

    print("Conversion completed successfully.")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Convert BAM file")
    parser.add_argument('bam_file', help='BAM file path')
    parser.add_argument('threads', type=int, help='Number of threads')
    parser.add_argument('output_prefix', help='Output prefix')
    parser.add_argument('reference_fasta', help='Reference genome FASTA file')
    parser.add_argument('input_sra', help='Input SRA number')
    args = parser.parse_args()

    convert_bam(args.bam_file, args.threads, args.output_prefix, args.reference_fasta, args.input_sra)