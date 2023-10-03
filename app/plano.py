import Levenshtein as lev
import blosum as bl
import math as ma


class SequenceProcessor:
    def __init__(
        self,
        kmerSize: int,
        sequence: str,
        blosumVersion: int = 62,
        logarithmBase: float = 2,
    ) -> None:
        self.kmerSize = kmerSize
        self.sequence = sequence
        self.blosumDict = dict(bl.BLOSUM(blosumVersion))
        self.minVariation = round(ma.log(self.kmerSize, logarithmBase))

    def __estimated_mutation_value(self, subsequence1, subsequence2) -> int:
        levDistance = lev.hamming(subsequence1, subsequence2)
        value = -10
        if levDistance <= self.minVariation:
            diffs = [
                self.blosumDict[subsequence1[i]][subsequence2[i]]
                for i in range(0, self.kmerSize)
                if subsequence1[i] != subsequence2[i]
            ]
            value = sum(diffs)
        return value

    def __create_graph(self, kmers, nodes) -> list:
        links = []
        for i in range(0, len(kmers) - 1):
            act, pos = kmers[i], kmers[i + 1]
            similarity = 0
            if act == pos:
                similarity = 0
            else:
                similarity = self.__estimated_mutation_value(act, pos)

            links.append([i, i + 1, similarity, nodes[act], nodes[pos]])
        return links

    def getRepeats(self) -> list | None:
        repeats = [] 
        for k in range(0, self.kmerSize):
            kmers = [
                self.sequence[i : i + self.kmerSize]
                for i in range(k, len(self.sequence), self.kmerSize)
                if len(self.sequence[i : i + self.kmerSize]) == self.kmerSize
            ]

            nodes = list(set(kmers))
            nodes = {nodes[j]: j for j in range(0, len(nodes))}
            links = self.__create_graph(kmers, nodes)
            nodes = {v: k for k, v in nodes.items()}
            repeats.extend(self.define_repeats(links, nodes)) 
        
        repeats 

    def define_repeats(self, links, nodes:dict):
        repeat = ""
        repeats = [] 
        ban  = False
        initialPosition = None 
        for i in range(len(links) - 1):
            current_sublist = links[i]
            next_sublist = links[i + 1]

            if current_sublist[1] == next_sublist[0] and current_sublist[2] >= 0:
                if not initialPosition: 
                    initialPosition = current_sublist[1]
                
                repeat = repeat + nodes.get(current_sublist[3])
                ban = True 
            elif ban: 
                repeat = repeat + nodes.get(current_sublist[3])
                repeats.append(repeat)
                initialPosition = None
                repeat = "" 
                ban  = False  

        return repeats 