import pandas as pd
import numpy as np 

class DeBruijnGraph:
  def __init__(self, kmer_size:int, sequence:str) -> None:
    self.kmer_size = kmer_size
    self.sequence = "".join(sequence.split())
    self.links = None
    self.kmers = []
    self.cycles = []

  def __build_graph(self, kmers:dict):
    columns = ['start', 'destine', 'link']
    link_number = 1
    links = []
    for k, v in kmers.items():
        pref, suf = v[:-1], v[1:]
        links.append([pref, suf, link_number])
        link_number += 1
    return pd.DataFrame(links, columns=columns)

  def __tarjan_find_cycles(self, links, current_link, max_link):

    tmp_cycles = {} 
    visited = {}
    cycles = []

    while current_link:
        start = current_link[0]
        destine = current_link[1]
        link_number = current_link[2]

        # Revisar la distancia máxima de separación
        if start in visited: # si ya se visitió el prefijo  
          if link_number - visited[start] <= 6: # determinar el tamañao del ciclo 
              if start not in tmp_cycles: # si no tiene creado un temporal de ciclos 
                tmp_cycles[start] = [[visited[start], link_number]] # se crea el temporal 
              else: 
                tmp_cycles[start].append([visited[start], link_number])   # si no, se extendiende la lista temporal de ciclos 
          else:
              if start in tmp_cycles: 
                if (len(tmp_cycles[start]) > 1 or tmp_cycles[start][0][0] + 1 == tmp_cycles[start][0][1]):  
                  cycles.extend(tmp_cycles[start]) # si el ciclo es demasiado grande, se guarda la lista de ciclos para ese prefijo 
                # tmp_cycles[start] = [[visited[start], link_number]] # se reinicia la temporal de ciclos - revisar si vacía o con el actual  
                del tmp_cycles[start]

        visited[start] = link_number

        if link_number < max_link:

            next_link = links[(links['start'] == destine) & (links['link'] == link_number + 1)].values.tolist()[0]
            current_link = next_link

        else:
            current_link = None
            break
       
    # si en los temporales quedan ciclos se evaluan y se guardan 
    for _, v in tmp_cycles.items(): 
      if len(v) == 1 and v[0] not in cycles and v[0][0] + 1 == v[0][1]: # ciclos de tipo NNN o KKK
        cycles.extend(v)
      if len(v) > 1:
        cycles_to_add = [ c for c in v if c not in cycles] 
        cycles.extend(cycles_to_add)    
      
    return cycles


  def __find_cycles(self, links: pd.DataFrame):
    # Encuentra el índice de la fila donde "Columna1" es mínimo
    indice_fila_minimo = links['link'].idxmin()

    # Usar el índice para obtener la fila completa
    fila_minimo = list(links.loc[indice_fila_minimo])

    # Encuentra el índice de la fila donde "Columna1" es máximo
    indice_fila_maximo = links['link'].idxmax()

    # Usar el índice para obtener la fila completa
    fila_maximo = list(links.loc[indice_fila_maximo])

    return self.__tarjan_find_cycles(links, fila_minimo, fila_maximo[2])

  def __get_cleanedRepeats(self, dirty_repeats):
    repeats = []
    merged_interval = None

    for interval in dirty_repeats:
        if not merged_interval:
            merged_interval = interval
        elif interval[0] <= merged_interval[1]:
            merged_interval[1] = max(merged_interval[1], interval[1])
        elif interval[0] == merged_interval[1] + 1: 
            merged_interval[1] = interval[1] 
        else:
            repeats.append(merged_interval)
            merged_interval = interval

    if merged_interval:
        repeats.append(merged_interval)

    return repeats

  def __establish_repeats(self, repeatStartEnd, sequence):
        repeatStrings = list()

        for repeat in repeatStartEnd:
            start = repeat[0] - 1
            end = repeat[1] + 1
            if end - start >= 5: 
              sequence = "".join(self.sequence[start:end])
              repeatStrings.append([sequence, start, end + 1])

        # Crear un DataFrame
        df = pd.DataFrame(repeatStrings, columns=["TR", "start", "end"])
        return df
  
  def __order_cycles(self, data): 
    cycles_ordered = [] 
    if len(data) > 1: 
      # Convertir la lista a un array de NumPy
      array_data = np.array(data)
      
      # Ordenar el array basado en el primer valor de cada sublista
      array_ordered_data = array_data[array_data[:,0].argsort()]
      
      # Convertir el array ordenado de nuevo a una lista
      cycles_ordered = array_ordered_data.tolist()

    return cycles_ordered

  def proccess(self):
    self.kmers = {i: self.sequence[i:i+self.kmer_size] for i in range(len(self.sequence)-self.kmer_size+1) if len(self.sequence[i:i+self.kmer_size]) == self.kmer_size}
    self.links = self.__build_graph(self.kmers)
    self.cycles = self.__find_cycles(self.links)
    self.cycles = self.__order_cycles(self.cycles)
    self.cycles = self.__get_cleanedRepeats(self.cycles)
    return self.__establish_repeats(self.cycles, self.sequence)