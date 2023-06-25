import subprocess

def execute_minimap2(reference_genome, fasta_file, threads):
    # Execute minimap2 command using subprocess
    minimap2_cmd = f"minimap2 -ax map-hifi -t {threads} --sam-hit-only --secondary=no {reference_genome} {fasta_file} > mapped.sam"
    try:
        subprocess.run(minimap2_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing minimap2 command: {e}")
