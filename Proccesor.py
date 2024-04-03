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
    parser.add_argument("--output_folder_path", type=validate_file_path, required=True,
                        help="Path to the output folder")
    parser.add_argument("--input_fasta_path", type=validate_file_path, required=True,
                        help="Path to the input fasta file")
    parser.add_argument("--k_mer_size", type=int, default=3,
                        help="Size of k-mer greater or equal to 3. (default: 3).")
    parser.add_argument("--sml", type=int, default=14,
                        help="SML (default: 14)")
    parser.add_argument("--minimal_size", type=int, default=9,
                        help="Minimal size of sequence found (possible tandem repeat) (default: 9)")
    parser.add_argument("--umbral", type=int, default=2,
                        help="Variation of between sml's of the same tandem repeat (default: 2)")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of workers for parallel processing (default: 4s)")

    args = parser.parse_args()

    # Validate additional parameters
    if args.k_mer_size <= 2:
        raise ValueError("k-mer size must be greater than three \" k-mer-size >=3 \"")
    if args.minimal_size < 9:
        raise ValueError("minimal_size size must be greater than 9 \" minimal_size >= 9 \"")
    if args.umbral < 0:
        raise ValueError("k-umbral size must be greater or equal than zero \" umbral >= 0 \"")
    if args.sml <= 0:
        raise ValueError("SML must be greater than zero \" sml > 0 \"")
    if args.workers <= 0:
        raise ValueError("Number of workers must be greater than zero")

    # Process sequences
    proccesor = FastaProccesor()
    proccesor.process_sequences_in_parallel(output_folder_path=args.output_folder_path,
                                                 input_fasta_path=args.input_fasta_path,
                                                 k=args.k_mer_size,
                                                 sml=args.sml,
                                                 minimal_size = args.minimal_size, 
                                                 umbral = args.umbral, 
                                                 workers=args.workers)
