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

TH_score_type = ["stable","raise","reduce"]

for s in T_vs_N_sample:
    os.system("mkdir -p Output_gene_log2Foldchange_TH_score/"+s+"-HiC")
    T_vs_N_combine_log2Foldchange = []
    T_vs_N_combine_log2Foldchange.append(["group","log2Foldchange"])
    for t in range(len(TH_score_type)):
        with open ("Input_TH_score_type_log2Foldchange/"+s+"-HiC/T_vs_N_"+TH_score_type[t]) as f:
            lines = f.readlines()
            for i in range(len(lines)):
                line = lines[i].strip()
                sep_line = line.split("\t")
                T_vs_N_combine_log2Foldchange.append([TH_score_type[t],sep_line[9]])
    savetxt("Output_gene_log2Foldchange_TH_score/"+s+"-HiC/T_vs_N_combine_log2Foldchange.csv",T_vs_N_combine_log2Foldchange)
    

for s in M_vs_T_sample:
    os.system("mkdir -p Output_gene_log2Foldchange_TH_score/"+s+"-HiC")
    M_vs_T_combine_log2Foldchange = []
    M_vs_T_combine_log2Foldchange.append(["group","log2Foldchange"])
    for t in range(len(TH_score_type)):
        with open ("Input_TH_score_type_log2Foldchange/"+s+"-HiC/M_vs_T_"+TH_score_type[t]) as f:
            lines = f.readlines()
            for i in range(len(lines)):
                line = lines[i].strip()
                sep_line = line.split("\t")
                M_vs_T_combine_log2Foldchange.append([TH_score_type[t],sep_line[9]])
    savetxt("Output_gene_log2Foldchange_TH_score/"+s+"-HiC/M_vs_T_combine_log2Foldchange.csv",M_vs_T_combine_log2Foldchange)
