import numpy as np
import os

def savetxt(filename, x):
    np.savetxt(filename, x, delimiter='\t', fmt='%s')

os.system("mkdir -p Output_plot_compartment_track")


os.system("pyGenomeTracks --tracks Input_compartment_tracks/tracks.ini --region chr20:0-63025520 -o Output_plot_compartment_track/chr20_compartment.svg --fontSize 10 --trackLabelHAlign left --width 30 --height 40")

