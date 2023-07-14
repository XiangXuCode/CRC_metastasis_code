import pycircos
import matplotlib.pyplot as plt
import collections
from matplotlib.backends.backend_pdf import PdfPages

Samples=["CRC-01-T-HiC", "CRC-02-M-HiC", "CRC-02-T-HiC", "CRC-03-LNM-HiC", "CRC-03-M-HiC", "CRC-03-T-HiC", "CRC-04-M-HiC", "CRC-04-T-HiC", "CRC-05-M-HiC", "CRC-05-T-HiC", "CRC-10-T-HiC", "CRC-11-T-HiC", "CRC-14-M-HiC", "CRC-14-T-HiC"]
colorlist = ["#ff8a80", "#ff80ab", "#ea80fc", "#b388ff", "#8c9eff", "#82b1ff", "#84ffff", "#a7ffeb", "#b9f6ca",
             "#ccff90", "#f4ff81", "#ffff8d", "#ffe57f", "#ffd180", "#ff9e80", "#bcaaa4", "#eeeeee", "#b0bec5",
             "#ff5252", "#ff4081", "#e040fb", "#7c4dff", "#536dfe", "#448aff", "#18ffff", "#64ffda", "#69f0ae",
             "#b2ff59", "#eeff41", "#ffff00", "#ffd740", "#ffab40", "#ff6e40", "#a1887f", "#e0e0e0", "#90a4ae"]

for sample in Samples:
    Garc    = pycircos.Garc
    Gcircle = pycircos.Gcircle

    #Set chromosomes
    circle = Gcircle(figsize=(10,10))
    with open("./3_circos/3_1_chromosome_general.csv") as f:
        f.readline()
        i = 0
        for line in f:
            line   = line.rstrip().split(",")
            name   = line[0]
            length = int(line[-1])
            arc    = Garc(arc_id=name, size=length, interspace=2,
                          raxis_range=(850,900), labelposition=80,
                          facecolor=colorlist[i],
                          label_visible=True)
            circle.add_garc(arc)
            i += 1
    circle.set_garcs(0,360)

    for arc_id in circle.garc_dict:
        circle.tickplot(arc_id, raxis_range=(900,915), tickinterval=20000000)

    values_all   = []
    arcdata_dict = collections.defaultdict(dict)
    with open("1_3_SV_all/%s.csv"%sample) as f:
        f.readline()
        for line in f:
            line  = line.rstrip().split(",")
            name1  = line[0]
            start1 = int(line[1])-1
            end1   = int(line[2])
            name2  = line[3]
            start2 = int(line[4])-1
            end2   = int(line[5])
            source = (name1, start1, end1, 850)
            destination = (name2, start2, end2, 850)
            circle.chord_plot(source, destination, facecolor=circle.garc_dict[name1].facecolor)
    pdf = PdfPages('./3_circos/3_3_all_%s.pdf'%sample)
    pdf.savefig()
    pdf.close()
    plt.close()
