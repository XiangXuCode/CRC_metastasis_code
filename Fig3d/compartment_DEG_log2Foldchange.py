import numpy as np
import os
import pandas as pd

def savetxt(filename, x):
    np.savetxt(filename, x, delimiter=',', fmt='%s')

chrname = []
for i in range(1,23):
    chrname.append("chr"+str(i))
chrname.append("chrX")

T_vs_N_sample = ["CRC-01"]
M_vs_T_sample = ["CRC-02"]

compartment_type = ["common_B","diff_A_to_B","diff_B_to_A","common_A"]

combine_significant_gene = []
with open("Input_DEseq2/significant_geneName") as f:
    lines = f.readlines()
    for i in range(len(lines)):
        line = lines[i].strip()
        combine_significant_gene.append(line.split('"')[1])
        
for s in T_vs_N_sample:
    os.system("mkdir -p Output_compartment_gene_log2Foldchange/"+s+"-HiC")
    T_vs_N_combine_log2Foldchange = []
    T_vs_N_combine_log2Foldchange.append(["group","log2Foldchange"])
    for t in range(len(compartment_type)):
        with open ("Input_compartment_DEG_log2Foldchange/"+s+"-HiC/T_vs_N_"+compartment_type[t]) as f:
            lines = f.readlines()
            for i in range(len(lines)):
                line = lines[i].strip()
                sep_line = line.split("\t")
                if sep_line[6] in combine_significant_gene:
                    T_vs_N_combine_log2Foldchange.append([compartment_type[t],sep_line[7]])
    savetxt("Output_compartment_gene_log2Foldchange/"+s+"-HiC/T_vs_N_combine_log2Foldchange.csv",T_vs_N_combine_log2Foldchange)
    

for s in M_vs_T_sample:
    os.system("mkdir -p Output_compartment_gene_log2Foldchange/"+s+"-HiC")
    M_vs_T_combine_log2Foldchange = []
    M_vs_T_combine_log2Foldchange.append(["group","log2Foldchange"])
    for t in range(len(compartment_type)):
        with open ("Input_compartment_DEG_log2Foldchange/"+s+"-HiC/M_vs_T_"+compartment_type[t]) as f:
            lines = f.readlines()
            for i in range(len(lines)):
                line = lines[i].strip()
                sep_line = line.split("\t")
                if sep_line[6] in combine_significant_gene:
                    M_vs_T_combine_log2Foldchange.append([compartment_type[t],sep_line[7]])
    savetxt("Output_compartment_gene_log2Foldchange/"+s+"-HiC/M_vs_T_combine_log2Foldchange.csv",M_vs_T_combine_log2Foldchange)
