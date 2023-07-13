#!/bin/bash
# bash ./unisort.sh -i <input.fasta> -o <output_prefix>
# Check if seqkit is present in PATH
if ! command -v seqkit &> /dev/null; then
    echo "seqkit not found, please install."
    exit 1
fi

while getopts ":i:o:" opt; do
  case $opt in
    i) input="$OPTARG"
    ;;
    o) output_prefix="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ -z "$input" ] || [ -z "$output_prefix" ]; then
  echo "Usage: bash unisort.sh -i <input.fasta> -o <output_prefix>"
  exit 1
fi

temp_file="${input}.int"
awk_script="/^>/{print \">${output_prefix}_\" ++i; next}{print}"

echo "seqkit found, proceeding to unisort"
seqkit sort -l -r "$input" > "$temp_file" &&
awk "$awk_script" "$temp_file" > "${output_prefix}.fasta" &&
rm "$temp_file"

