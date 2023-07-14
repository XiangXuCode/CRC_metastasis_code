import numpy as np
import os
import matplotlib.pyplot as plt

def savetxt(filename, x):
    np.savetxt(filename, x, delimiter='\t', fmt='%s')

chrname = []
for i in range(1,23):
    chrname.append("chr"+str(i))
chrname.append("chrX")

T_vs_N_sample = ["CRC-03"]
M_vs_T_sample = ["CRC-03"]
M_vs_N_sample = ["CRC-03"]

for s in T_vs_N_sample:
    os.system("mkdir -p Output_compartment_change/"+s+"-HiC")
    N = []
    with open("Input_compartment/"+s+"-N-HiC/compartment.txt") as f:
        lines = f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            N.append([sep_line[2],sep_line[3],sep_line[4],sep_line[7]])
    T = []
    with open("Input_compartment/"+s+"-T-HiC/compartment.txt") as f:
        lines = f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            T.append([sep_line[2],sep_line[3],sep_line[4],sep_line[7]])    

    T_vs_N = []
    T_vs_N_common = []
    T_vs_N_common_A = []
    T_vs_N_common_B = []
    T_vs_N_diff = []
    T_vs_N_diff_A_to_B = []
    T_vs_N_diff_B_to_A = []
    for j in range(len(lines)-1):
        if N[j][3] != "NA" and T[j][3] != "NA":
            T_vs_N.append([N[j][0],N[j][1],N[j][2],str(float(T[j][3])-float(N[j][3]))])
            if float(N[j][3])>0 and float(T[j][3])>0:
                T_vs_N_common_A.append([N[j][0],N[j][1],N[j][2]])
                T_vs_N_common.append([N[j][0],N[j][1],N[j][2]])
            elif float(N[j][3])<0 and float(T[j][3])<0:
                T_vs_N_common_B.append([N[j][0],N[j][1],N[j][2]])
                T_vs_N_common.append([N[j][0],N[j][1],N[j][2]])
            elif float(N[j][3])>0 and float(T[j][3])<0:
                T_vs_N_diff_A_to_B.append([N[j][0],N[j][1],N[j][2]])
                T_vs_N_diff.append([N[j][0],N[j][1],N[j][2]])
            elif float(N[j][3])<0 and float(T[j][3])>0:
                T_vs_N_diff_B_to_A.append([N[j][0],N[j][1],N[j][2]])
                T_vs_N_diff.append([N[j][0],N[j][1],N[j][2]])  
    
    data = [len(T_vs_N_diff_B_to_A),len(T_vs_N_common_A),len(T_vs_N_common_B),len(T_vs_N_diff_A_to_B)]
    labels = ["B_to_A","Common_A","Common_B","A_to_B"]
    explode = [0.3,0.02,0.02,0.3]
    colors = ["#c21f30","#eea2a4","#baccd9","#15559a"]
    plt.pie(data,labels=labels,colors = colors,autopct = '%1.2f%%',explode = explode)
    plt.title(s+"-HiC")
    plt.savefig("Output_compartment_change/"+s+"-HiC/T_vs_N_pie.svg")
    plt.cla()
    

for s in M_vs_T_sample:
    os.system("mkdir -p Output_compartment_change/"+s+"-HiC")
    T = []
    with open("Input_compartment/"+s+"-T-HiC/compartment.txt") as f:
        lines = f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            T.append([sep_line[2],sep_line[3],sep_line[4],sep_line[7]])
    M = []
    with open("Input_compartment/"+s+"-M-HiC/compartment.txt") as f:
        lines = f.readlines()
        for i in range(1,len(lines)):
            line = lines[i].strip()
            sep_line = line.split("\t")
            M.append([sep_line[2],sep_line[3],sep_line[4],sep_line[7]])    

    M_vs_T = []
    M_vs_T_common = []
    M_vs_T_common_A = []
    M_vs_T_common_B = []
    M_vs_T_diff = []
    M_vs_T_diff_A_to_B = []
    M_vs_T_diff_B_to_A = []
    for j in range(len(lines)-1):
        if T[j][3] != "NA" and M[j][3] != "NA":
            M_vs_T.append([T[j][0],T[j][1],T[j][2],str(float(M[j][3])-float(T[j][3]))])
            if float(T[j][3])>0 and float(M[j][3])>0:
                M_vs_T_common_A.append([T[j][0],T[j][1],T[j][2]])
                M_vs_T_common.append([T[j][0],T[j][1],T[j][2]])
            elif float(T[j][3])<0 and float(M[j][3])<0:
                M_vs_T_common_B.append([T[j][0],T[j][1],T[j][2]])
                M_vs_T_common.append([T[j][0],T[j][1],T[j][2]])
            elif float(T[j][3])>0 and float(M[j][3])<0:
                M_vs_T_diff_A_to_B.append([T[j][0],T[j][1],T[j][2]])
                M_vs_T_diff.append([T[j][0],T[j][1],T[j][2]])
            elif float(T[j][3])<0 and float(M[j][3])>0:
                M_vs_T_diff_B_to_A.append([T[j][0],T[j][1],T[j][2]])
                M_vs_T_diff.append([T[j][0],T[j][1],T[j][2]])


    data = [len(M_vs_T_diff_B_to_A),len(M_vs_T_common_A),len(M_vs_T_common_B),len(M_vs_T_diff_A_to_B)]
    labels = ["B_to_A","Common_A","Common_B","A_to_B"]
    explode = [0.3,0.02,0.02,0.3]
    colors = ["#c21f30","#eea2a4","#baccd9","#15559a"]
    plt.pie(data,labels=labels,colors = colors,autopct = '%1.2f%%',explode = explode)
    plt.title(s+"-HiC")
    plt.savefig("Output_compartment_change/"+s+"-HiC/M_vs_T_pie.svg")
    plt.cla()

