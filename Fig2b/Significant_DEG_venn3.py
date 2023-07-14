import os
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

def savetxt(filename,x):
    np.savetxt(filename,x,delimiter = '\t',fmt='%s')

T_vs_N = []
T_vs_N_gene = []
f = open("Input_DEseq2/T_vs_N/adjqvalue_foldchange.csv","r")
lines = f.readlines()
T_vs_N.append(lines[0].strip().split(","))
for i in range(1,len(lines)):
    line = lines[i].strip()
    sep_line = line.split(",")
    if sep_line[6] != "NA" and sep_line[2] != "NA":
        if float(sep_line[6]) <= 0.05 and abs(float(sep_line[2])) >= 1:
            T_vs_N.append(sep_line)
            T_vs_N_gene.append(sep_line[0])
f.close()


M_vs_T = []
M_vs_T_gene = []
f = open("Input_DEseq2/M_vs_T/adjqvalue_foldchange.csv","r")
lines = f.readlines()
M_vs_T.append(lines[0].strip().split(","))
for i in range(1,len(lines)):
    line = lines[i].strip()
    sep_line = line.split(",")
    if sep_line[6] != "NA" and sep_line[2] != "NA":
        if float(sep_line[6]) <= 0.05 and abs(float(sep_line[2])) >= 1:
            M_vs_T.append(sep_line)
            M_vs_T_gene.append(sep_line[0])
f.close()


M_vs_N = []
M_vs_N_gene = []
f = open("Input_DEseq2/M_vs_N/adjqvalue_foldchange.csv","r")
lines = f.readlines()
M_vs_N.append(lines[0].strip().split(","))
for i in range(1,len(lines)):
    line = lines[i].strip()
    sep_line = line.split(",")
    if sep_line[6] != "NA" and sep_line[2] != "NA":
        if float(sep_line[6]) <= 0.05 and abs(float(sep_line[2])) >= 1:
            M_vs_N.append(sep_line)
            M_vs_N_gene.append(sep_line[0])
f.close()


overlap = list((set(T_vs_N_gene).intersection(set(M_vs_T_gene))).intersection(set(M_vs_N_gene)))

combine = list((set(T_vs_N_gene).union(set(M_vs_T_gene))).union(set(M_vs_N_gene)))

os.system("mkdir -p Output_venn")

venn3(subsets = [set(T_vs_N_gene),set(M_vs_T_gene),set(M_vs_N_gene)],set_labels = ("T_vs_N","M_vs_T","M_vs_N"),set_colors = ("#ff0000","#ffbb00","#00ffef"))
plt.savefig("Output_venn/venn3.svg")

