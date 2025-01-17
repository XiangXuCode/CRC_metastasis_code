import os

chrname = []
for i in range(1,23):
    chrname.append("chr"+str(i))
chrname.append("chrX")

os.system("mkdir -p Output_hicPlotTADs")
    
os.system("hicPlotTADs --tracks Input_TAD_tracks/tracks.ini -o Output_hicPlotTADs/hic_track.svg --region chr8:100000000-140000000")    



