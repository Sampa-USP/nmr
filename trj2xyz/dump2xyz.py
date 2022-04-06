import tqdm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sea
import pandas as pd 

class dump:
    def __init__(self, time_step, xyz_df):
        self.time_step = time_step
        self.len = len(xyz_df)
        self.xyz_df = xyz_df

    def get_xyz_type(self, type_num):
        df_selected = self.xyz_df[self.xyz_df[1] == type_num]
        return df_selected

def read_dumps(file_path,frames, extra_line_num,verbose=False):
    line_pointer = 0
    num_of_atom_df = pd.read_csv(file_path, skiprows=line_pointer + 3, \
                            nrows=1, header=None)
    num_of_atom = int(num_of_atom_df[0])
    dump_list = []
    
    aux = 0
    
    for _ in tqdm.tqdm(range(frames)):
        try:
            time_step_df = pd.read_csv(file_path, skiprows=line_pointer + 1, \
                                    nrows=1, header=None)
        except:
            break
            
        time_step = int(time_step_df[0])

        xyz_df = pd.read_csv(file_path, skiprows=line_pointer + 9, \
                                   nrows=num_of_atom, sep=' ', header=None)
        xyz_df = xyz_df.iloc[:,:]
        line_pointer += (extra_line_num + num_of_atom)
        dump_list.append(dump(time_step, xyz_df))

    return dump_list 


def save_all_dumps(list_of_dumps,path_to_save):
    aux = 0
    for dump in tqdm.tqdm(list_of_dumps):
        dump.xyz_df[2] = dump.xyz_df[2].map({2:"H",1:"O"})
        dump.xyz_df = teste[0].xyz_df.sort_values(by=[0])

        array_xyz = np.array(dump.xyz_df)[:,2:7]

        texto = str(len(array_xyz))+"\n"+"\n"

        for linha in array_xyz:
            texto += ' '.join(list(linha.astype("str"))) + "\n"
            #texto+ =aux
        texto = texto[:-1]

        file = open(path_to_save+str(aux).zfill(len(str(int(len(list_of_dumps)))))+".xyz",'w')
        file.write(texto)
        file.close()
        
        aux+=1
#!ls *.GraphGeod
#!pwd
