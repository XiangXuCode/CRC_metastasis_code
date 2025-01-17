import numpy as np
import os

def savetxt(filename, x):
    np.savetxt(filename, x, delimiter='\t', fmt='%s')

chrname = []
for i in range(1,23):
    chrname.append("chr"+str(i))
chrname.append("chrX")

without_metastasis_T_vs_N = ["CRC-01-HiC","CRC-10-HiC","CRC-11-HiC"]
metastasis_T_vs_N = ["CRC-02-HiC","CRC-03-HiC","CRC-04-HiC","CRC-05-HiC","CRC-14-HiC"]
metastasis_M_vs_T = ["CRC-02-HiC","CRC-03-HiC","CRC-04-HiC","CRC-05-HiC","CRC-14-HiC"]
metastasis_M_vs_N = ["CRC-02-HiC","CRC-03-HiC","CRC-04-HiC","CRC-05-HiC","CRC-14-HiC"]

os.system("mkdir -p Output_compartment_stability")

compartment_stability = []
compartment_stability.append(["group","SOD"])

for s in without_metastasis_T_vs_N:
    with open("Input_compartment_change/"+s+"/T_vs_N_common") as f:
        lines = f.readlines()
        common = len(lines)
    with open("Input_compartment_change/"+s+"/T_vs_N_diff") as f:
        lines = f.readlines()
        diff = len(lines)
    compartment_stability.append(["without_metastasis_T_vs_N",str(common/(common+diff))])

for s in metastasis_T_vs_N:
    with open("Input_compartment_change/"+s+"/T_vs_N_common") as f:
        lines = f.readlines()
        common = len(lines)
    with open("Input_compartment_change/"+s+"/T_vs_N_diff") as f:
        lines = f.readlines()
        diff = len(lines)
    compartment_stability.append(["metastasis_T_vs_N",str(common/(common+diff))])

for s in metastasis_M_vs_T:
    with open("Input_compartment_change/"+s+"/M_vs_T_common") as f:
        lines = f.readlines()
        common = len(lines)
    with open("Input_compartment_change/"+s+"/M_vs_T_diff") as f:
        lines = f.readlines()
        diff = len(lines)
    compartment_stability.append(["metastasis_M_vs_T",str(common/(common+diff))])

for s in metastasis_M_vs_N:
    with open("Input_compartment_change/"+s+"/M_vs_N_common") as f:
        lines = f.readlines()
        common = len(lines)
    with open("Input_compartment_change/"+s+"/M_vs_N_diff") as f:
        lines = f.readlines()
        diff = len(lines)
    compartment_stability.append(["metastasis_M_vs_N",str(common/(common+diff))])

savetxt("Output_compartment_stability/compartment_stability.txt",compartment_stability)

