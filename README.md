# DeBruijn-Based Tandem Repeats Finder

The DeBruijn-Based Tandem Repeats Finder is a bioinformatics tool designed to identify tandem repeats within protein sequences. Unlike DNA-based tandem repeat finders, which search for repeated DNA sequences, this tool focuses on identifying patterns of repeated amino acids within proteins. It utilizes the DeBruijn graph algorithm, a method commonly used in sequence analysis, to efficiently detect tandem repeats in protein sequences. This tool is valuable for studying the structural and functional implications of tandem repeats in proteins, which can play important roles in protein folding, stability, and function.

### Dockerized app  
docker run  -v /resources/:/var/resources/ bruijn_proccesor python /app/Proccesor.py --input_fasta_path ./var/resources/uniprotkb_taxonomy_id_137071_2023_09_13.fasta --output_folder_path ./var/resources