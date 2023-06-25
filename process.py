import os
import subprocess

class SRAImportResult:
    def __init__(self, fastq_file_path, fasta_file_path):
        self.fastq_file_path = fastq_file_path
        self.fasta_file_path = fasta_file_path

def process_sra(sra_number):
    # Execute prefetch command using subprocess
    prefetch_cmd = f"prefetch {sra_number}"
    try:
        subprocess.run(prefetch_cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error executing prefetch command: {e}")
        return

    # Extract the directory name
    directory_name = sra_number

    # Check if the directory exists
    if os.path.exists(directory_name):
        # Change working directory to the prefetch-created directory
        os.chdir(directory_name)

        # Execute fasterq-dump command using subprocess
        fasterq_cmd = f"fasterq-dump {sra_number}"
        try:
            subprocess.run(fasterq_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing fasterq-dump command: {e}")
            return

        # Capture the generated .fastq file path
        fastq_file = f"{sra_number}.fastq"
        fastq_file_path = os.path.join(os.getcwd(), fastq_file)

        # Execute seqtk command to convert .fastq to .fasta
        fasta_file = f"{sra_number}.fasta"
        fasta_file_path = os.path.join(os.getcwd(), fasta_file)
        seqtk_cmd = f"seqtk seq -A {fastq_file_path} > {fasta_file_path}"
        try:
            subprocess.run(seqtk_cmd, shell=True, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error executing seqtk command: {e}")
            return

        # Change back to the original working directory
        os.chdir('..')

        # Return the SRA import result
        return SRAImportResult(fastq_file_path, fasta_file_path)
    else:
        print(f"Directory '{directory_name}' does not exist.")
        return
