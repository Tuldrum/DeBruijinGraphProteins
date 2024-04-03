import os
import time
import pandas as pd
import re 
import concurrent.futures

from datetime import datetime
from Bio import SeqIO, SeqRecord

from app.DeBruijnGraph import DeBruijnGraph
from app.GraphProccesor import GraphProccesor


class FastaProccesor(object):

    def __get_outputfile_name(self, local_file_fasta, output_folder_path):
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        pattern = r"[\\\/]([^\\\/]+)\.\w+$"
        match = re.search(pattern, local_file_fasta)
        file_name = 'results'
        if match: 
            file_name = match.group(1)
        output_file_path = os.path.join(output_folder_path, f'{file_name}_{timestamp}.csv') 
        return output_file_path

    def __get_seq_name(self, name: str):
        return name.split("|")[-1]

    def process_sequence(self, record: SeqRecord, k: int = 3, sml: int = 9, minimal_size: int = 9, umbral: int = 2): 
        return self.__process_sequence(record, k, sml, minimal_size, umbral)

    def __process_sequence(self, record: SeqRecord, k: int = 3, sml: int = 9, minimal_size: int = 9, umbral: int = 2):
        seq = str(record.seq)
        start = time.time()
        bruijn = DeBruijnGraph(sequence=seq, k=k)
        proccesor = GraphProccesor(
            deBruijinGraph=bruijn, sml=sml, minimal_size = minimal_size, umbral=umbral, 
        )
        repeats = proccesor.proccess()
        end = time.time()
        duration = end - start
        # repeats["SequenceId"] = self.__get_seq_name(record.id)
        repeats["SequenceId"] = record.id
        repeats["Duration"] = duration
        print(repeats.shape[0], record.id)
        return repeats

    def __validate_files(self, input_fasta_path, output_folder_path):
        ban = True
        try:
            if input_fasta_path and output_folder_path:
                # Check if input FASTA file exists
                if not os.path.exists(input_fasta_path):
                    print("Input FASTA file does not exist.")
                    ban = False

                # Create output directory if it doesn't exist
                output_directory = os.path.dirname(output_folder_path)
                if not os.path.exists(output_directory) and ban:
                    os.makedirs(output_directory)
                    print(f"Output directory '{output_directory}' created.")

            elif not input_fasta_path:
                ban = False
                print("Please provide a fasta file path")
            elif not output_folder_path:
                ban = False
                print("Please provide a path to store results")

        except Exception as e:
            ban = False
            print(e)

        return ban

    def process_sequences_in_parallel(self, input_fasta_path, output_folder_path, k=3, sml=9, workers=4, minimal_size = 9, umbral = 14):
        if self.__validate_files(input_fasta_path=input_fasta_path, output_folder_path=output_folder_path):
            all_repeats = None
            with concurrent.futures.ProcessPoolExecutor(
                max_workers=workers
            ) as executor:
                records = list(SeqIO.parse(input_fasta_path, "fasta"))
                results = executor.map(
                    self.process_sequence,
                    records,
                    [k]*len(records), 
                    [sml]*len(records), 
                    [minimal_size]*len(records), 
                    [umbral]*len(records)
                )

                for repeats in results:
                    if all_repeats is not None:
                        all_repeats = pd.concat([all_repeats, repeats], axis=0)
                    else:
                        all_repeats = repeats

            if all_repeats is not None:
                file_path = self.__get_outputfile_name(local_file_fasta=input_fasta_path, output_folder_path=output_folder_path)
                all_repeats.to_csv(file_path, index=False)
                print(f"Results saved to '{output_folder_path}'.")
