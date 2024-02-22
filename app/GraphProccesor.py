from DeBruijnGraph import DeBruijnGraph
from iteration_utilities import flatten
import pandas as pd
import numpy as np


class GraphProccesor:

    def __init__(self, deBruijinGraph: DeBruijnGraph, minimal_residues: int) -> None:
        self.__deBruijinGraph = deBruijinGraph
        self.__minimal_residues = minimal_residues

    def __establish_cycles(self):

        links = self.__deBruijinGraph.graph
        # Obtener la tercera columna
        tercera_columna = links[:, 2]

        # Encontrar el mayor y menor elemento de la tercera columna
        max_link = np.max(tercera_columna)
        current_link = np.min(tercera_columna)

        # Encontrar los índices donde los valores son iguales al máximo y al mínimo
        current_link_indices = np.where(tercera_columna == current_link)[0]
        max_link_indices = np.where(tercera_columna == max_link)[0]

        # Obtener las filas correspondientes a los índices encontrados
        current_link = links[current_link_indices][0]
        max_link = links[max_link_indices][0]

        tmp_cycles = {}
        visited = {}
        cycles = []

        while current_link is not None:
            start = current_link[0]
            destine = current_link[1]
            link_number = current_link[2]

            # Revisar la distancia máxima de separación
            if start in visited:
                if link_number - visited[start] <= self.__minimal_residues:
                    if start not in tmp_cycles:
                        tmp_cycles[start] = [[visited[start], link_number]]
                    else:
                        tmp_cycles[start].append([visited[start], link_number])
                else:
                    if start in tmp_cycles:
                        if (
                            len(tmp_cycles[start]) > 1
                            or tmp_cycles[start][0][0] + 1 == tmp_cycles[start][0][1]
                        ):
                            cycles.extend(tmp_cycles[start])
                        del tmp_cycles[start]

            visited[start] = link_number

            if link_number < max_link[2]:
                next_link = links[
                    (links[:, 0] == destine) & (links[:, 2] == link_number + 1)
                ][0]
                current_link = next_link
            else:
                current_link = None
                break

        # Si en los temporales quedan ciclos, se evalúan y se guardan
        for _, v in tmp_cycles.items():
            if len(v) == 1 and v[0] not in cycles and v[0][0] + 1 == v[0][1]:
                cycles.extend(v)
            if len(v) > 1:
                cycles_to_add = [c for c in v if c not in cycles]
                cycles.extend(cycles_to_add)

        return cycles

    def __joining_cycles(self, cycles):
        cycles = sorted(cycles, key=lambda x: x[0])
        previous_cycle: list = None
        cycles_found = []
        for cycle in cycles:
            cycle = [cycle]
            if previous_cycle is None:
                previous_cycle = cycle
            else:
                a_0, a_1 = np.min(previous_cycle), np.max(previous_cycle)
                b_0, _ = np.min(cycle), np.max(cycle)

                if b_0 > a_1:

                    if len(previous_cycle) > 1:
                        cycles_found.append(previous_cycle)
                    previous_cycle = cycle

                elif a_0 <= b_0 <= a_1:
                    previous_cycle.extend(cycle)

        if previous_cycle:
            cycles_found.append(previous_cycle)

        return cycles_found

    def __get_str_repeats(self, cycles) -> list:
        str_cycles: list = []
        k_left = self.__deBruijinGraph.kmer_size - 2
        for cycle in cycles:
            flattened_cylce = list(flatten(cycle))
            x0, x1 = min(flattened_cylce), max(flattened_cylce) + k_left
            if x1 - x0 > 3: 
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
        return cycles
