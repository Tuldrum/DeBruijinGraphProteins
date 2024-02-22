import json
import numpy as np


class DeBruijnGraph:
    """
    Class to represent a De Bruijn graph.
    """

    def __init__(self, sequence: str = "", k: int = 9) -> None:
        """
        Initialize the De Bruijn graph with a given sequence and k-mer length.

        Args:
            sequence (str): The input DNA sequence.
            k (int): The length of k-mers to consider.
        """
        self.__sequence = sequence
        self.__kmer_size = k
        self.__normalized_nodes = {}
        self.__graph = self.__generate_graph()

    @property
    def graph(self):
        return self.__graph.copy()

    @property
    def kmer_size(self):
        return self.__kmer_size

    @kmer_size.setter
    def kmer_size(self, k: int):
        self.__kmer_size = k
        self.__graph = self.__generate_graph()

    @property
    def sequence(self):
        return self.__sequence

    @sequence.setter
    def sequence(self, seq: str):
        self.__sequence = seq
        self.__graph = self.__generate_graph()

    def __str__(self) -> str:
        return json.dumps(self.__graph)

    def __generate_graph(self):
        """
        Generate the De Bruijn graph based on the sequence and k-mer length.

        Returns:
            links: The De Bruijn graph representation.
        """
        sequence = self.sequence
        k = self.kmer_size
        cnt = 0
        links = list()
        for i in range(len(self.sequence) - (k - 1)):
            pref, suf = sequence[i : k + i - 1], sequence[i + 1 : k + i]
            if pref not in self.__normalized_nodes:
                cnt += 1
                self.__normalized_nodes[pref] = cnt
            if suf not in self.__normalized_nodes:
                cnt += 1
                self.__normalized_nodes[suf] = cnt

            pref, suf = self.__normalized_nodes[pref], self.__normalized_nodes[suf]
            links.append([pref, suf, i])
        return np.array(links)

    def get_sub_sequence(self, x0, x1):
        subsequence = "Error"
        if x1 > x0 and self.__sequence and len(self.sequence) > x1 + 1:
            subsequence = self.__sequence[x0 : x1 + 1]

        return subsequence
