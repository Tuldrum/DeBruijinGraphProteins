from Bio.Align import PairwiseAligner

# Define las secuencias
seq1 = "SEQ1"
seq2 = "SEQ2"

# Crea un objeto PairwiseAligner
aligner = PairwiseAligner()

# Define la matriz de sustitución, por ejemplo, BLOSUM62
aligner.substitution_matrix = "BLOSUM62"

# Realiza el alineamiento
alignments = aligner.align(seq1, seq2)

# Imprime el mejor alineamiento
best_alignment = alignments[0]
print("Alineamiento:")
print(best_alignment)