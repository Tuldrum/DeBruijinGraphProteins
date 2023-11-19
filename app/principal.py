from Bio import SeqIO
import time
from FASTAhandler import FastaDownloader 
from BruijnGraph import DeBruijnGraph
import pandas as pd 
from datetime import datetime
import concurrent.futures
import re

local_file_fasta = ".\\app\\resources\\downloaded.fasta"

def getname(name): 
    ## patron = r'[^|]+\|[^|]+\|(.*)'
    ## resultado = re.search(patron, name)
    return name.split("|")[-1]

def process_sequence(record):
    seq = str(record.seq)
    # if record.id.startswith('tr'): 
    start = time.time()
    bruijn = DeBruijnGraph(sequence=seq, kmer_size=3)
    print(record.id)
    repeats = bruijn.proccess()
    end = time.time()
    duration = end - start 
    repeats['SequenceId'] = getname(record.id)
    repeats['Duration'] = duration 
    print(repeats.size, record.id)
    return repeats

if __name__ == "__main__":
    all_repeats = None

    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        records = list(SeqIO.parse(local_file_fasta, "fasta"))
        results = executor.map(process_sequence, records)

        for repeats in results:
            if all_repeats is not None: 
                all_repeats = pd.concat([all_repeats, repeats], axis=0)
            else: 
                all_repeats = repeats 

    file_name = '.\\app\\resources\\datos_' + datetime.now().strftime('%Y%m%d_%H%M%S') + '.csv'

    if all_repeats is not None: 
        all_repeats.to_csv(file_name, index=False)
