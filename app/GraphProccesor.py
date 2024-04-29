from iteration_utilities import flatten
import pandas as pd
import numpy as np
import copy
from app.DeBruijnGraph import DeBruijnGraph


class GraphProccesor:

    def __init__(self, 
                deBruijinGraph: DeBruijnGraph, 
                sml: int, minimal_size: int, umbral:int, 
                max_dist_pos:int, min_cycles_joined:int) -> None:
        self.__deBruijinGraph = deBruijinGraph
        self.__sml = sml
        self.__minimal_size = minimal_size
        self.__umbral = umbral
        self.__repeats_found = None 
        self._k = self.__deBruijinGraph.kmer_size
        self.__max_dist_pos = max_dist_pos 
        self.__min_cycles_joined = min_cycles_joined

    @property
    def deBruijnGraph(self) -> DeBruijnGraph : 
        return copy.deepcopy(self.__deBruijinGraph) 

    @deBruijnGraph.setter
    def deBruijnGraph(self, newGraph: DeBruijnGraph): 
        if isinstance(newGraph, DeBruijnGraph): 
            self.__deBruijinGraph = newGraph
        else: 
            raise TypeError("Is not a DeBruijnGraph instance class")

    def __get_max_and_min_links(self):
        links = self.__deBruijinGraph.graph
        # Obtener la tercera columna
        tercera_columna = links[:, 2]

        # Encontrar el mayor y menor elemento de la tercera columna
        max_link = np.max(tercera_columna)
        # Encontrar el elemento más pequeño de la tercera columna
        current_link = np.min(tercera_columna)

        # Encontrar los índices donde los valores son iguales al máximo y al mínimo
        current_link_indices = np.where(tercera_columna == current_link)[0]
        max_link_indices = np.where(tercera_columna == max_link)[0]

        # Obtener las filas correspondientes a los índices encontrados
        current_link = links[current_link_indices][0]
        max_link = links[max_link_indices][0]

        return current_link, max_link

    def __validate_umbral(self, tmp_dist_cycles, new_dist_cycle): 
        min_dist = abs(min(tmp_dist_cycles) - new_dist_cycle)
        max_dist = abs(max(tmp_dist_cycles) - new_dist_cycle)
        return True if min_dist <= self.__umbral and max_dist <= self.__umbral else False
    

    def __check_and_update_cycles(self, start, link_number, visited, tmp_cycles, cycles, tpm_dist_cycles):
        ban = True 
        if (link_number - 1)  - (visited[start] + self._k-2) <= self.__sml:
            if start not in tmp_cycles:
                tmp_cycles[start] = [[visited[start], link_number]]
                tpm_dist_cycles[start] = [abs(visited[start] - link_number)]
            else:
                dist = abs(visited[start] - link_number)
                if self.__validate_umbral(tpm_dist_cycles[start], dist): 
                    tmp_cycles[start].append([visited[start], link_number])
                    tpm_dist_cycles[start].append(dist)
                else: 
                    ban = False 
        else: 
            ban = False

        if ban is False: 
            if start in tmp_cycles:
                if len(tmp_cycles[start]) > 1 or tmp_cycles[start][0][0] + 1 == tmp_cycles[start][0][1]:
                    cycles.extend(tmp_cycles[start])
                del tmp_cycles[start]
                del tpm_dist_cycles[start]

    def __get_next_link(self, links, destine, link_number, max_link):
        next_link = None
        if link_number < max_link[2]:
            next_link = links[(links[:, 0] == destine) & (links[:, 2] == link_number + 1)][0]
        return next_link

    def __add_remaining_cycles(self, tmp_cycles, cycles):
        for _, v in tmp_cycles.items():
            # Si solo hay un subciclo en la lista y no está en los ciclos encontrados y es de tipo NNN 
            # if len(v) == 1 and v[0] not in cycles and v[0][0] + 1 == v[0][1]: 
            #     print(v)
            #     cycles.extend(v)
            if len(v) > 1: # [[], []]
                cycles_to_add = [c for c in v if c not in cycles]
                cycles.extend(cycles_to_add)

    def __establish_cycles(self):
        links = self.__deBruijinGraph.graph
        current_link, max_link = self.__get_max_and_min_links() 

        tmp_cycles = {}
        tpm_dist_cycles = {}
        visited = {}
        cycles = []

        while current_link is not None:
            start = current_link[0]
            destine = current_link[1]
            link_number = current_link[2]

            # Revisar la distancia máxima de separación
            if start in visited:
                self.__check_and_update_cycles(start, link_number, visited, tmp_cycles, cycles, tpm_dist_cycles)

            visited[start] = link_number
            current_link = self.__get_next_link(links, destine, link_number, max_link)

        # Si en los temporales quedan ciclos, se evalúan y se guardan
        self.__add_remaining_cycles(tmp_cycles, cycles)
        return cycles

    def __joining_cycles(self, cycles):
        cycles = sorted(cycles, key=lambda x: x[0])
        previous_cycle: list | None = None
        cycles_found = []
        cycles_joined = 0 
        last_position_vs: int | None = None

        for cycle in cycles: 
            if previous_cycle is None: 
                previous_cycle = [np.min(cycle), np.max(cycle) + (self._k-2)]
                last_position_vs = previous_cycle[0]
            else: 
                
                a_0, a_1 = np.min(previous_cycle), np.max(previous_cycle) # ciclo ant 
                b_0, b_1 = np.min(cycle), np.max(cycle) + (self._k-2) # ciclo nuevo 

                if b_0 > a_1: 
                    if cycles_joined >= self.__min_cycles_joined: 
                        cycles_found.append(previous_cycle)
                    previous_cycle = [b_0, b_1]
                    cycles_joined = 0

                elif a_0 <= b_0 <= a_1 and b_0 - last_position_vs <= self.__max_dist_pos:
                    previous_cycle = [a_0, b_1 if b_1 > a_1 else a_1]
                    cycles_joined+=1 

                last_position_vs = b_0
                
        if previous_cycle and cycles_joined >= self.__min_cycles_joined:
            print(previous_cycle)
            cycles_found.append(previous_cycle)

        return cycles_found

    def __joining_cycles_2(self, cycles):
        cycles = sorted(cycles, key=lambda x: x[0])
        previous_cycle: list = None
        cycles_found = []
        k_const = (self._k-2)
        counter_cycles = 1
        for cycle in cycles:
            if previous_cycle is None:
                previous_cycle = [cycle[0], cycle[1] + k_const]
            else:
                a_0, a_1 = np.min(previous_cycle), np.max(previous_cycle) # ciclo ant 
                b_0, b_1 = np.min(cycle), np.max(cycle) + (self._k-2) # ciclo nuevo 

                if b_0 > a_1: 
                    if counter_cycles > 1: 
                        cycles_found.append(previous_cycle)
                    previous_cycle = [b_0, b_1]
                    counter_cycles = 1 

                elif a_0 <= b_0 <= a_1:
                    previous_cycle = [a_0, b_1 if b_1 > a_1 else a_1]
                    counter_cycles = counter_cycles + 1

        if previous_cycle and counter_cycles > 1:
            cycles_found.append(previous_cycle)

        return cycles_found

    def __get_str_repeats(self, cycles) -> list:
        str_cycles: list = []
        for cycle in cycles:
            x0, x1 = min(cycle), max(cycle) 
            if x1 - x0 + 1 >= self.__minimal_size:
                str_cycle = self.__deBruijinGraph.get_sub_sequence(x0=x0, x1=x1)
                str_cycles.append([str_cycle, x0, x1])
        return str_cycles

    def __to_dataframe(self, cycles):
        df = pd.DataFrame(cycles, columns=["TR", "start", "end"])
        return df

    def proccess(self):
        cycles = self.__establish_cycles()
        cycles = self.__joining_cycles(cycles)
        cycles = self.__get_str_repeats(cycles)
        cycles = self.__to_dataframe(cycles)
        self.__repeats_found = cycles
        return self.__repeats_found

    @property
    def repeats_found(self):
        return self.__repeats_found  