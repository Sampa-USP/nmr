import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large'}


pylab.rcParams.update(params)

import sys
import seaborn as sea
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.ticker import PercentFormatter
import glob
import numpy as np
import pandas as pd
import tqdm
import argparse

sys.path.append("/home/jvilela/python_libs/")
import chemnetworks
import lammps

def load_graph_path(path_chem):
    L = []
    for files in glob.glob(path_chem):
        for files2 in glob.glob(files+"/*.GraphGeod"):
            L+=[files2]
    return L

def load_xyz(path_chem,interval=10):
    L = []
    for files in glob.glob(path_chem):
        for files2 in glob.glob(files+"/*.xyz")[::interval]:
            L+=[files2]
    return L

def load_dicts(paths):
    full_dict = {}
    full_array = {}
    count = 0
    for path in paths:
        lh_ = chemnetworks.read_GraphGeod(path)
        dict_ = lh_.bond_dict
        full_dict[str(count)] = dict_
        
        df = pd.DataFrame(lh_.matriz,columns="id_a_1 id_a_2 iX iY iZ tp_a_1 tp_a_2 d_AA angulo".split(" "))
        
        full_array[str(count)] = df

        count+=1
    return full_dict,full_array



#parser = argparse.ArgumentParser(description='Description of your program')
#parser.add_argument('-f','--ff_folder', help='get the main name of ff folders', required=True)
#print(args = vars(parser.parse_args()))
ff_main_name = "tip3p"



all_paths = load_graph_path(f"../{ff_main_name}/CHEM_out/*")
full_dicts_2005,full_array = load_dicts(all_paths)
all_xyz_paths = load_xyz(f"../{ff_main_name}/tjr4xyz")

count_1 = 0
count_2 = 0
count_3 = 0
count_4 = 0
count_5 = 0

counter = 0

keys = np.array(list(full_array.keys())).astype(str)

for xyz_files in tqdm.tqdm(all_xyz_paths):
    teste = open(xyz_files,"r").read()
    s_array = teste.split("\n")[2:]
    id_dict = {}
    for i in range(len(s_array[::3])):
        id_dict[i] = s_array[3*i:3*i+3]
        
    partial_array = full_array[keys[counter]]

    for atom1 in np.unique(partial_array['id_a_1']):
        lhs = np.sum(partial_array['id_a_1'].values==atom1)
        
        hbons = full_array[keys[counter]]["id_a_2"][full_array[keys[counter]]["id_a_1"]==atom1].values
        
        if lhs==1:
            count_1 +=1
            file = open(f"lh_clusters_xyz/{ff_main_name}/1/{count_1}.xyz","w")
            file.write(f"{3*(lhs+1)}\n\n")
            file.write("\n".join(id_dict[int(atom1)-1]))
            for h2os in hbons:
                file.write("\n")
                file.write("\n".join(id_dict[int(h2os)-1]))
            file.close()

        if lhs==2:
            count_2 +=1
            file = open(f"lh_clusters_xyz/{ff_main_name}/2/{count_2}.xyz","w")
            file.write(f"{3*(lhs+1)}\n\n")
            file.write("\n".join(id_dict[int(atom1)-1]))
            for h2os in hbons:
                file.write("\n")
                file.write("\n".join(id_dict[int(h2os)-1]))
            file.close()
            
            
        if lhs==3:
            count_3 +=1
            file = open(f"lh_clusters_xyz/{ff_main_name}/3/{count_3}.xyz","w")
            file.write(f"{3*(lhs+1)}\n\n")
            file.write("\n".join(id_dict[int(atom1)-1]))
            for h2os in hbons:
                file.write("\n")
                file.write("\n".join(id_dict[int(h2os)-1]))
            file.close()
        
        if lhs==4:
            count_4 +=1
            file = open(f"lh_clusters_xyz/{ff_main_name}/4/{count_4}.xyz","w")
            file.write(f"{3*(lhs+1)}\n\n")
            file.write("\n".join(id_dict[int(atom1)-1]))
            for h2os in hbons:
                file.write("\n")
                file.write("\n".join(id_dict[int(h2os)-1]))
                file.write("\n")
            file.close()
            
            
        if lhs==5:
            count_5 +=1
            file = open(f"lh_clusters_xyz/{ff_main_name}/5/{count_5}.xyz","w")
            file.write(f"{3*(lhs+1)}\n\n")
            file.write("\n".join(id_dict[int(atom1)-1]))
            for h2os in hbons:
                file.write("\n")
                file.write("\n".join(id_dict[int(h2os)-1]))
            file.close()
    counter+=1
    
out_put_thermo = {}

for number_of_h_bonds in range(1,6):
    all_xyzs = glob.glob(f"lh_clusters_xyz/{ff_main_name}/{number_of_h_bonds}/*")
    
    print(f"lh_clusters_xyz/{ff_main_name}/{number_of_h_bonds}/*")
    out_put_thermo[str(number_of_h_bonds)] = []

    charge0 = -1.1128 #tip4p/2005 O charge
    charge1 = 0.5564 #tip4p/2005 H charge

    for xyz_s in range(len(all_xyzs)):
        #print(xyz_s)
        #xyz_s=0
        file = open(all_xyzs[xyz_s],"r").read().split("\n")


        topol_init = f"""LAMMPS Description

            {str(int(int(file[0])))}  atoms
            {str(int(int(file[0])*2/3))}  bonds
            {str(int(int(file[0])*1/3))}  angles
            0  dihedrals
            0  impropers

            2  atom types
            1  bond types
            1  angle types

          0.0 43.7400 xlo xhi
          0.0 43.7400 ylo yhi
          0.0 43.7400 zlo zhi

        Masses

        1 15.9994  # O
        2 1.008  # H

        Atoms  # full

        """

        mol_count = 0
        end_topol = ""
        final_bonds = []
        end_bonds = ""

        final_angs = []
        end_angs = ""
        
        file = [x for x in file if x]
        file = [file[0]]+[""]+file[1:]
        
        for elem in range(len(file[2::3])):
            splited_string_0 = file[elem*3+2].split(" ")
            splited_string_H_1 = file[elem*3+3].split(" ")
            splited_string_H_2 = file[elem*3+4].split(" ")

            string_0 = [str(elem*3+1)] + [str(elem+1)] + [str(1)] + [str(charge0)] + (" ".join(splited_string_0[-3:])).split(" ")
            string_h1 = [str(elem*3+2)] + [str(elem+1)] + [str(2)] + [str(charge1)] + (" ".join(splited_string_H_1[-3:])).split(" ")    
            string_h2 = [str(elem*3+3)] + [str(elem+1)] + [str(2)] + [str(charge1)] + (" ".join(splited_string_H_2[-3:])).split(" ") 

            end_topol += (' '.join(string_0) + '\n')
            end_topol += (' '.join(string_h1) + '\n')
            end_topol += (' '.join(string_h2) + '\n')

            end_bonds += ' '.join([str(elem*2+1)] + [str(1)] + [str(elem*3+1)] + [str(elem*3+2)])+"\n"
            end_bonds += ' '.join([str(elem*2+2)] + [str(1)] + [str(elem*3+1)] + [str(elem*3+3)])+"\n"

            end_angs += ' '.join([str(1+elem)] + [str(1)] + [str(elem*3+2)] + [str(elem*3+1)] + [str(elem*3+3)])+"\n"

        end_topol = topol_init+end_topol+"\n"+'Bonds\n\n'+end_bonds+'\n'+'Angles\n\n'+end_angs+'\n'

        topo_file = open(f"topo_file_{ff_main_name}_tmp.data","w")
        topo_file.write(end_topol)
        topo_file.close()

        str_system_in_init = f"""
        # -- Default styles (for solo "SPCE" water) --
        units        real
        atom_style   full
        pair_style   lj/cut/coul/long 12.0
        pair_modify  tail yes mix arithmetic
        bond_style   harmonic
        angle_style  harmonic
        kspace_style pppm 0.00001
        #pair_modify  mix arithmetic  # LEAVE THIS UNSPECIFIED!

        read_data topo_file_{ff_main_name}_tmp.data

        bond_coeff   1         450.0   0.9572
        angle_coeff  1       55.0    104.52
        pair_coeff   1 1  0.102   3.188  
        pair_coeff   2 2  0.0     0.0
        pair_coeff 1 2 0.0 0.0
        group spce type  1  2

        group central_h20 molecule 1
        group others_molecule molecule {' '.join(np.array(range(2,2+1)).astype(str))}

        thermo_style custom evdwl ecoul elong epair pe
        thermo 1

        run 1
        """
        lmp = lammps.lammps()
        lmp.commands_list(str_system_in_init.split("\n"))
        out_put_thermo[str(number_of_h_bonds)] += [[lmp.get_thermo("evdwl")],[lmp.get_thermo("ecoul")],[lmp.get_thermo("elong")],[lmp.get_thermo("epair")],[lmp.get_thermo("pe")]]
        lmp.close()

    ar = np.array(out_put_thermo[str(number_of_h_bonds)])
    ar = ar.reshape((len(ar)//5,5))
    np.savetxt(f"energy_arrays/{ff_main_name}/ar_{number_of_h_bonds}.txt",ar)
    
    #if count==