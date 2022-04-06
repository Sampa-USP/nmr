import numpy as np


class read_GraphGeod():
    def __init__(self,file_path):
        """id átomo 1 | id átomo 2 | iX | iY | iZ | timo átomo 1 | tipo átomo 2 | distância (em angstrongs) | ângulo ($\theta$)"""
        self.file = file_path
        self.reader()
        self.plot_values()
        
        
    def reader(self):

        r_file = open(self.file,"r").read().split("\n")
        matriz = [r_file[i].split(" ") for i in range(len(r_file)-1)]
        matriz = np.array(matriz).astype(float)

        matriz2 = matriz.copy()
        aux = matriz[:,0]
        matriz2[:,0] = matriz[:,1]
        matriz2[:,1] = aux

        Matriz = np.concatenate((matriz,matriz2))

        matriz = Matriz.copy()

        self.matriz = matriz[np.argwhere(matriz[:,0]!=matriz[:,1]).flatten()]

        graph = {}



        for atom1 in np.unique(matriz[:,0]):
            graph[str(atom1)] = len(np.where(matriz[:,0]==atom1)[0].flatten())

        self.graph = graph

    def plot_values(self):
        self.bond_dict = {'0':np.sum(np.array(list(self.graph.values()))==0),
             '1':np.sum(np.array(list(self.graph.values()))==1),
             '2':np.sum(np.array(list(self.graph.values()))==2),
             '3':np.sum(np.array(list(self.graph.values()))==3),
             '4':np.sum(np.array(list(self.graph.values()))==4),
             '5':np.sum(np.array(list(self.graph.values()))==5),
             '6':np.sum(np.array(list(self.graph.values()))==6)}
