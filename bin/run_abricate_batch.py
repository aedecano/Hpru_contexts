import os
import subprocess
import argparse

def run_abricate(abricate_path, fasta_directory, output_directory):
    # Create output directory if it doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Iterate through each FASTA file in the directory
    for fasta_file in os.listdir(fasta_directory):
        if fasta_file.endswith(".fasta") or fasta_file.endswith(".fa"):
            # Construct full paths
            fasta_path = os.path.join(fasta_directory, fasta_file)
            output_path = os.path.join(output_directory, f"{os.path.splitext(fasta_file)[0]}_abricate_results.txt")

            # Run Abricate using subprocess
            command = f"{abricate_path} {fasta_path} > {output_path}"
            
            # Adjust the command based on Abricate's actual command-line parameters
            # The example assumes the input FASTA file is passed directly, adjust as needed

            subprocess.run(command, shell=True)

    print("Abricate analysis completed.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Abricate on FASTA files in a directory.")
    parser.add_argument("--abricate-path", required=True, help="Path to the Abricate script.")
    parser.add_argument("--fasta-directory", required=True, help="Directory containing FASTA files.")
    parser.add_argument("--output-directory", required=True, help="Output directory for Abricate results.")
    
    args = parser.parse_args()

    run_abricate(args.abricate_path, args.fasta_directory, args.output_directory)
