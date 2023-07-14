import numpy as np
import os
import pandas as pd

def savetxt(filename, x):
    np.savetxt(filename, x, delimiter='\t', fmt='%s')

chrname = []
for i in range(1,23):
    chrname.append("chr"+str(i))
chrname.append("chrX")

file_path = "Input_dense_matrix/"
dirs = os.listdir(file_path)

resolution = 40000

threshold = 0.6

min_square_size = 6*resolution
max_square_size = 6000000

width = [1,2,3]

print("parameters: resolution="+str(resolution)+" threshold="+str(threshold)+" min_square_size="+str(min_square_size)+" max_square_size="+str(max_square_size)+" width="+str(width))

for s in dirs:
    print(s)
    os.system("mkdir -p Output_square_domain/"+str(resolution)+"/"+s)
    square_domain = []
    square_domain_juicer = []
    square_anchor = []
    all_ratio_list = []
    all_ratio_list.append(["c1","start","end","c2","start","end","min_ratio","inside_x_positive_ratio","inside_y_positive_ratio","outside_x_positive_ratio","outside_y_positive_ratio"])
    for c in chrname:
        print(c)
        dense_matrix = []
        with open("Input_dense_matrix/"+s+"/"+str(resolution)+"/"+c+".matrix") as f:
            lines = f.readlines()
            for i in range(len(lines)):
                line = lines[i].strip()
                sep_line = line.split("\t")
                dense_matrix.append([float(x) for x in sep_line])
        h = np.array(dense_matrix)

        for a in range(len(lines)):
            for b in range(int(min_square_size/resolution),int(max_square_size/resolution)+1):
                gap = 0
                if a>=1 and a+b+1 <= len(lines):
                    gap = 1
                    for e in width:
                        gap = gap*float(dense_matrix[a+e][a+e])*float(dense_matrix[a+b-1-e][a+b-1-e])*float(dense_matrix[a-1][a-1])*float(dense_matrix[a+b][a+b])
                if a>=1 and a+b+1 <= len(lines) and gap!=0.0 and a+b+max(width) <= len(lines):
                    diff_width_ratio_min = []
                    diff_width_ratio = []
                    for g in width:
                        inside_x_positive = 0
                        inside_x_pairs = b+1-2*g
                        inside_y_positive = 0
                        inside_y_pairs = b+1-2*g

                        if np.mean(h[a:a+g,a:a+g])>dense_matrix[a+g][a+g]:
                            inside_x_positive = inside_x_positive +0.5
                        if np.mean(h[a+b-g:a+b,a+b-g:a+b])>dense_matrix[a+b-1-g][a+b-1-g]:
                            inside_y_positive = inside_y_positive +0.5
                        if np.mean(h[a:a+g,a+b-g:a+b])>dense_matrix[a+g][a+b-1-g]:
                            inside_x_positive = inside_x_positive +0.5
                            inside_y_positive = inside_y_positive +0.5
                        for d in range(b-2*g):  
                            if np.mean(h[a:a+g,a+g+d])>dense_matrix[a+g][a+d+g]:
                                inside_x_positive = inside_x_positive +1            
                            if np.mean(h[a+d+g,a+b-g:a+b])>dense_matrix[a+d+g][a+b-1-g]:
                                inside_y_positive = inside_y_positive +1 

                        outside_x_positive = 0
                        outside_x_pairs = b+1
                        outside_y_positive = 0
                        outside_y_pairs = b+1
                        if dense_matrix[a-1][a-1]<np.mean(h[a:a+g,a:a+g]):
                            outside_x_positive = outside_x_positive +0.5
                        if dense_matrix[a+b][a+b]<np.mean(h[a+b-g:a+b,a+b-g:a+b]):
                            outside_y_positive = outside_y_positive +0.5
                        if dense_matrix[a-1][a+b]<np.mean(h[a:a+g,a+b-g:a+b]):
                            outside_x_positive = outside_x_positive +0.5
                            outside_y_positive = outside_y_positive +0.5
                        for d in range(b):  
                            if dense_matrix[a-1][a+d]<np.mean(h[a:a+g,a+g+d]):
                                outside_x_positive = outside_x_positive +1            
                            if dense_matrix[a+d][a+b]<np.mean(h[a+d+g,a+b-g:a+b]):
                                outside_y_positive = outside_y_positive +1 

                        diff_width_ratio_min.append(min(float(inside_x_positive/inside_x_pairs),float(inside_y_positive/inside_y_pairs),float(outside_x_positive/outside_x_pairs),float(outside_y_positive/outside_y_pairs)))
                        diff_width_ratio.append([float(inside_x_positive/inside_x_pairs),float(inside_y_positive/inside_y_pairs),float(outside_x_positive/outside_x_pairs),float(outside_y_positive/outside_y_pairs)])


                    best_width_num = diff_width_ratio_min.index(max(diff_width_ratio_min))
                    best_width = width[best_width_num]

                    all_ratio_list.append([c,str(a*resolution),str((a+best_width)*resolution),c,str((a+b-best_width)*resolution),str((a+b)*resolution),max(diff_width_ratio_min),diff_width_ratio[best_width_num][0],diff_width_ratio[best_width_num][1],diff_width_ratio[best_width_num][2],diff_width_ratio[best_width_num][3]])

                    if max(diff_width_ratio_min)>=threshold:
                            square_domain.append([c,str(a*resolution),str((a+b)*resolution)])
                            square_domain_juicer.append([c,str(a*resolution),str((a+b)*resolution),c,str(a*resolution),str((a+b)*resolution)])
                            square_anchor.append([c,str(a*resolution),str((a+best_width)*resolution),c,str((a+b-best_width)*resolution),str((a+b)*resolution)])

    savetxt("Output_square_domain/"+str(resolution)+"/"+s+"/square_domain.bed",square_domain)
    savetxt("Output_square_domain/"+str(resolution)+"/"+s+"/square_domain_juicer",square_domain_juicer)
    savetxt("Output_square_domain/"+str(resolution)+"/"+s+"/square_anchor",square_anchor)
    savetxt("Output_square_domain/"+str(resolution)+"/"+s+"/all_ratio_list",all_ratio_list)
    print("square_domain_founded = "+str(len(square_domain)))
                               
             
        



