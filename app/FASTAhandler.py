import requests
from Bio import SeqIO
import os 

class FastaDownloader: 

    def __init__(self, fasta_url:str, 
                 local_file_path:str = ".\\app\\resources\\downloaded.fasta") -> None:
        self.fasta_url = fasta_url 
        self.local_file_path = local_file_path
    
    def __download_fasta_file(self): 
        downloaded: bool = False 
        try:
            response = requests.get(self.fasta_url)

            # Check if the request was successful (status code 200)
            if response.status_code == 200:
                # Open the local file for writing in binary mode
                with open(self.local_file_path, "wb") as fasta_file:
                    # Write the content of the HTTP response to the local file
                    fasta_file.write(response.content)
                downloaded = True 
                print(f"FASTA file downloaded and saved to {self.local_file_path}")
            else:
                print(f"Failed to download FASTA file. HTTP status code: {response.status_code}")
        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")
        return downloaded 

    def get_fasta_records(self):
        if self.__download_fasta_file():
            if not os.path.isfile(self.local_file_path):
                raise FileNotFoundError(f"The downloaded file {self.local_file_path} does not exist.")
            
            return list(SeqIO.parse(self.local_file_path, "fasta"))
        else:
            raise Exception("No se pudo descargar el archivo")