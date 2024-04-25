import re
import statistics
import numpy as np 

def prefix_repetitions_distribution(s, k):
    repetitions = {}
    # pattern = rf'(?=(([a-zA-Z]{{{k}}}).*\2))'
    pattern = rf'(?=(([a-zA-Z]{{{k}}}).*))'
    matches = re.finditer(pattern, s)
    
    for match in matches:
        prefix = match.group(2)
        repetition_start = match.start()
        repetitions.setdefault(prefix, []).append(repetition_start)
    
    results = {}
    for prefix, positions in repetitions.items():
        if len(positions) >= 2: 
            distances = [(positions[i+1] - 1) - (positions[i] + (k-1)) for i in range(len(positions)-1)]
            # avg_distance = statistics.(distances) if distances else 0
            distancia_mas_frecuente = max(distances, key=distances.count)
            std_dev = statistics.stdev(distances) if len(distances) > 1 else 0
            mad = np.median(np.abs(distances - np.median(distances))) if len(distances) > 1 else 0
            num_repetitions = len(positions)
            results[prefix] = {'distancia_mas_frecuente': distancia_mas_frecuente, 'std_dev': std_dev, 'mad': mad, 'num_repetitions': num_repetitions}
    
    return results

# Ejemplo de uso
cadena = "LSFESKLNQGNKDILSLNEMNESELFLSFNDFYVYPSKFFLLSKKCKKIILLFLPKNVQSYEERLFL"
k = 2
resultado = prefix_repetitions_distribution(cadena, k)
for prefix, info in resultado.items():
    print(f"Prefijo: {prefix}, distancia_mas_frecuente: {info['distancia_mas_frecuente']}, Desviación estándar: {info['std_dev']}, MAD: {info['mad']}, Cantidad de repeticiones: {info['num_repetitions']}")