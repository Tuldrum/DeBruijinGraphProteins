import argparse
from app.FastaProcessor import FastaProccesor

def validate_file_path(path):
    # Basic validation for file path
    if not path:
        raise ValueError("File path cannot be empty")
    # Additional validation can be added here
    return path

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="DeBruijn-Based Tandem Repeats Finder")
    parser.add_argument("--output_csv_path", type=validate_file_path, required=True,
                        help="Path to the output CSV file")
    parser.add_argument("--input_fasta_path", type=validate_file_path, required=True,
                        help="Path to the input fasta file")
    parser.add_argument("--k_mer_size", type=int, default=3,
                        help="Size of k-mer (default: 3)")
    parser.add_argument("--maximum_residues", type=int, default=9,
                        help="Maximum number of residues (default: 9)")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of workers for parallel processing (default: 4s)")

    args = parser.parse_args()

    # Validate additional parameters
    if args.k_mer_size <= 2:
        raise ValueError("k-mer size must be greater than three \" k-mer-size >=3 \"")
    if args.maximum_residues <= 0:
        raise ValueError("Maximum residues must be greater than zero")
    if args.workers <= 0:
        raise ValueError("Number of workers must be greater than zero")

    # Process sequences
    proccesor = FastaProccesor()
    proccesor.process_sequences_in_parallel(output_csv_path=args.output_csv_path,
                                                 input_fasta_path=args.input_fasta_path,
                                                 k=args.k_mer_size,
                                                 maximum_residues=args.maximum_residues,
                                                 workers=args.workers)