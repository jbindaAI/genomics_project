import os
import zipfile
import pickle
import subprocess
import argparse
from Bio import SeqIO


def fetch_proteomes_ncbi_datasets(accession_file, output_dir):
    """
    Fetch proteomes from NCBI using the Datasets CLI for given Assembly Accessions.

    Parameters:
    - accession_file: str, Path to a text file containing Assembly Accessions.
    - output_dir: str, Directory to save the downloaded proteomes.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(accession_file, "r") as f:
        for line in f:
            # Parse Assembly Accession and optional description
            parts = line.strip().split(";")
            accession = parts[0].strip()
            description = parts[1].strip() if len(parts) > 1 else "Unknown organism"

            print(f"Fetching proteome for accession: {accession} ({description})")

            try:
                # Construct a command
                output_file = os.path.join(output_dir, f"{accession}.zip")
                command = [
                    "datasets", "download", "genome",
                    "accession", accession,
                    "--include", "protein",
                    "--filename", output_file
                ]
                
                # Execute the command
                subprocess.run(command, check=True)
                print(f"Saved proteome archive to {output_file}")

            except subprocess.CalledProcessError as e:
                print(f"Error fetching proteome for {accession}: {e}")


def extract_protein_faa(zip_dir, output_dir):
    """
    Extract the protein.faa file from each ZIP archive in the specified directory.

    Parameters:
    - zip_dir: str, Path to the directory containing ZIP archives.
    - output_dir: str, Path to the directory to save extracted protein.faa files.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    for filename in os.listdir(zip_dir):
        if filename.endswith(".zip"):
            zip_path = os.path.join(zip_dir, filename)
            try:
                with zipfile.ZipFile(zip_path, 'r') as z:
                    # Find the desired file in the archive
                    target_file = None
                    for file in z.namelist():
                        if file.endswith("/protein.faa"):
                            target_file = file
                            break

                    if target_file:
                        print(f"Extracting {target_file} from {filename}")
                        output_file_path = os.path.join(output_dir, os.path.basename(os.path.dirname(target_file)) + ".faa")
                        with z.open(target_file) as source, open(output_file_path, 'wb') as target:
                            target.write(source.read())
                        print(f"Saved to {output_file_path}")
                    else:
                        print(f"No protein.faa file found in {filename}")

            except zipfile.BadZipFile:
                print(f"Error: {filename} is not a valid zip file.")


def prepare_genome_names_map(accession_file, output_dir):
    """
    Create a mapping from genome IDs to natural names based on accesion file.
    
    Parameters:
    - accession_file: str, File with genome IDs and corresponding names.
    - output_dir: str, Directory to save a genome name map.
    - output_name: str, Name of the prepared genome name map.
    """

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    print("Preparing mapping from genome IDs to natural names.")

    name_map = {}
    with open(accession_file, "r") as f:
        for line in f:
            parts = line.strip().split(";")
            accession = parts[0].strip()
            name = parts[1].strip()
            name_map[accession] = name
    
    #with open(os.path.join(output_dir, "genomeID2name.pkl"), "wb") as f:
        #pickle.dump(name_map, f)
    return name_map


def prepare_genome_map(protein_files_dir, genomeID2name, output_dir, output_name="genome_map.pkl"):
    """
    Create a mapping from protein sequence IDs to genome IDs, genome names, and sequences based on protein FASTA files.
    
    Parameters:
    - protein_files_dir: str, Directory containing protein FASTA files.
    - genomeID2name: dict, Mapping of genome IDs to genome names.
    - output_dir: str, Directory to save the genome map.
    - output_name: str, Name of the file to save the genome map.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    print("Preparing mapping from protein IDs to genome IDs, genome names, and sequences.")

    genome_map = {}
    for file in os.listdir(protein_files_dir):
        if file.endswith(".faa"):
            genome_id = os.path.splitext(file)[0]
            genome_name = genomeID2name.get(genome_id, "Unknown Genome")

            file_path = os.path.join(protein_files_dir, file)
            for record in SeqIO.parse(file_path, "fasta"):
                prot_id = record.id  # Protein ID
                sequence = str(record.seq)  # Protein sequence
                genome_map[prot_id] = (genome_id, genome_name, sequence)

    output_path = os.path.join(output_dir, output_name)
    with open(output_path, "wb") as f:
        pickle.dump(genome_map, f)


def combine_fasta_files(input_dir:str, output_dir:str):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    combined_fasta = os.path.join(output_dir, "combined_proteins.faa")

    # Combine all FASTA files into one
    with open(combined_fasta, "w") as outfile:
        for file in os.listdir(input_dir):
            if file.endswith(".faa"):
                with open(os.path.join(input_dir, file), "r") as infile:
                    outfile.write(infile.read())


def is_done(path2check:str)->bool:
    if os.path.exists(path2check) and len(os.listdir(path2check)) > 0:
        return True
    else:
        return False
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script for fetching and processing proteomes.")
    parser.add_argument("--accession_filename", required=True, help="Name of an accession file.")
    args = parser.parse_args()

    BASENAME = args.accession_filename.split(".")[0]
    ACCESSION_FILEPATH = os.path.join("data_preparation/accession_ids", args.accession_filename)
    ARCHIVE_DIR = os.path.join("data_preparation/data/proteome_archives", BASENAME)
    FASTA_DIR = os.path.join("data_preparation/data/proteome_fasta", BASENAME)
    COMBINED_FASTA_DIR = os.path.join("data_preparation/data/combined_fasta", BASENAME)
    MAPS_DIR = os.path.join("data_preparation/data/maps", BASENAME)

    if not is_done(ARCHIVE_DIR):
        fetch_proteomes_ncbi_datasets(ACCESSION_FILEPATH, ARCHIVE_DIR)
    else:
        print("Data already downloaded!")
    
    if not is_done(FASTA_DIR):
        extract_protein_faa(ARCHIVE_DIR, FASTA_DIR)
    else:
        print("Data already preprocessed!")

    if not is_done(COMBINED_FASTA_DIR):
        combine_fasta_files(FASTA_DIR, COMBINED_FASTA_DIR)
    else:
        print("Fasta files already merged!")

    if not is_done(MAPS_DIR):
        genome_name_map = prepare_genome_names_map(ACCESSION_FILEPATH, MAPS_DIR)
        prepare_genome_map(FASTA_DIR, genome_name_map, MAPS_DIR, output_name="genome_map.pkl")
    else:
        print("Maps are already prepared!")
