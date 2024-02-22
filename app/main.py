from Bio import SeqIO
import time
from DeBruijnGraph import DeBruijnGraph
import pandas as pd 
from datetime import datetime
import concurrent.futures
from GraphProccesor import GraphProccesor 

local_file_fasta = "app\\resources\\uniprotkb_taxonomy_id_36329_2023_10_31.fasta"

def getFilename(filename:str): 
    return filename.replace('app\\resources\\', "").replace(".fasta", "")

def getname(name): 
    ## patron = r'[^|]+\|[^|]+\|(.*)'
    ## resultado = re.search(patron, name)
    return name.split("|")[-1]

def process_sequence(record):
    seq = str(record.seq)
    # if record.id.startswith('tr'): 
    start = time.time()
    bruijn = DeBruijnGraph(sequence=seq, k=3)
    # bruijn = DeBruijnGraph(sequence=record, k=3)
    proccesor = GraphProccesor(deBruijinGraph=bruijn, minimal_residues=9)
    # print(record.id)
    repeats = proccesor.proccess() 
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

    file_name = '.\\app\\resources\\datos_' + getFilename(local_file_fasta) + '_' + datetime.now().strftime('%Y%m%d_%H%M%S') + '.csv'

    if all_repeats is not None: 
        all_repeats.to_csv(file_name, index=False)