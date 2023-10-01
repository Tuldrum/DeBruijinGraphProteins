import requests

# Replace this URL with the URL of the FASTA file you want to download
fasta_url = "https://example.com/path/to/your.fasta"

# Specify the local file path where you want to save the downloaded FASTA file
local_file_path = "downloaded.fasta"

try:
    response = requests.get(fasta_url)

    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Open the local file for writing in binary mode
        with open(local_file_path, "wb") as fasta_file:
            # Write the content of the HTTP response to the local file
            fasta_file.write(response.content)
        print(f"FASTA file downloaded and saved to {local_file_path}")
    else:
        print(f"Failed to download FASTA file. HTTP status code: {response.status_code}")
except requests.exceptions.RequestException as e:
    print(f"An error occurred: {e}")




"""

 http://artemisa.unicauca.edu.co/~eumartinez/Proteinas/uniprotkb_taxonomy_id_137071_2023_09_13.fasta


from Bio import SeqIO

for index,record in enumerate(SeqIO.parse("uniprotkb_taxonomy_id_137071_2023_09_13.fasta", "fasta")):
    print(record.id)
    print(repr(record.seq))
    print(len(record))

"""