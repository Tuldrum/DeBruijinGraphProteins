from plano import SequenceProcessor 
from FASTAhandler import FastaDownloader 

if __name__ == '__main__':

    fasta_url = "http://artemisa.unicauca.edu.co/~eumartinez/Proteinas/uniprotkb_taxonomy_id_137071_2023_09_13.fasta" 
    local_file_path = "ruta_local_del_archivo.fasta"
    fasta_downloader = FastaDownloader(fasta_url)
    records = fasta_downloader.get_fasta_records()

    record = records[0]
    print(record)
    processor = SequenceProcessor(kmerSize=4, sequence=str(record.seq))
    processor.getRepeats() 
