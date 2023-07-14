import sys
from neoloop.visualize.core import *
import cooler

def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
mkdir('./2_1_heatmap')

def plot_heat(mcool, fragment):
    clr = cooler.Cooler('./0_0_mcool/%s.mcool::resolutions/50000' % mcool)
    translocation = 'C26     translocation,' + fragment
    vis = Triangle(clr, translocation, n_rows=6, figsize=(7, 5.3), track_partition=[5, 0.4, 0.8, 0.3, 0.3, 0.5],
                   correct='weight', space=0.03)
    vis.matrix_plot(vmin=0, vmax=0.001, cbr_fontsize=9, no_colorbar=True)
    vis.plot_chromosome_bounds(linewidth=2)
    vis.plot_genes(filter_=["ZMIZ1"], fontsize=10, release=75)
    mkdir('./2_1_heatmap/%s' % "ZMIZ1")
    vis.outfig('./2_1_heatmap/%s/%s_%s.pdf' % ("ZMIZ1", mcool, "ZMIZ1"), dpi=300)

Samples = ["CRC-02-N-HiC", "CRC-02-T-HiC", "CRC-02-M-HiC"]
for sample in Samples:
    plot_heat(sample, "3,26900000,+,10,70500000,- 3,4600000 10,92500000")

def plot_heat(mcool, fragment):
    clr = cooler.Cooler('./0_0_mcool/%s.mcool::resolutions/50000' % mcool)
    translocation = 'C26     translocation,' + fragment
    vis = Triangle(clr, translocation, n_rows=6, figsize=(7, 5.3), track_partition=[5, 0.4, 0.8, 0.3, 0.3, 0.5],
                   correct='weight', space=0.03)
    vis.matrix_plot(vmin=0, vmax=0.001, cbr_fontsize=9, no_colorbar=True)
    vis.plot_chromosome_bounds(linewidth=2)
    vis.plot_genes(filter_=["TTC28"], fontsize=10, release=75)
    mkdir('./2_1_heatmap/%s' % "TTC28")
    vis.outfig('./2_1_heatmap/%s/%s_%s.pdf' % ("TTC28", mcool, "TTC28"), dpi=300)

Samples = ["CRC-02-N-HiC", "CRC-02-T-HiC", "CRC-02-M-HiC"]
for sample in Samples:
    plot_heat(sample, "5,107200000,+,22,17900000,- 5,83500000 22,42700000")


# def plot_heat (mcool, fragment, gene):
def plot_heat(mcool, fragment):
    clr = cooler.Cooler('./0_0_mcool/%s.mcool::resolutions/50000' % mcool)
    translocation = 'C26     translocation,' + fragment
    vis = Triangle(clr, translocation, n_rows=6, figsize=(7, 5.3), track_partition=[5, 0.4, 0.8, 0.3, 0.3, 0.5],
                   correct='weight', space=0.03)
    vis.matrix_plot(vmin=0, vmax=0.001, cbr_fontsize=9, no_colorbar=True)
    vis.plot_chromosome_bounds(linewidth=2)
    vis.plot_genes(filter_=["CLN3"], fontsize=10, release=75)
    mkdir('./2_1_heatmap/%s' % "CLN3")
    vis.outfig('./2_1_heatmap/%s/%s_%s.pdf' % ("CLN3", mcool, "CLN3"), dpi=300)

Samples = ["CRC-05-N-HiC", "CRC-05-T-HiC", "CRC-05-M-HiC"]
for sample in Samples:
    plot_heat(sample, "16,31900000,+,22,32200000,- 16,23000000 22,36200000")


def plot_heat(mcool, fragment):
    clr = cooler.Cooler('./0_0_mcool/%s.mcool::resolutions/50000' % mcool)
    translocation = 'C26     translocation,' + fragment
    vis = Triangle(clr, translocation, n_rows=6, figsize=(7, 5.3), track_partition=[5, 0.4, 0.8, 0.3, 0.3, 0.5],
                   correct='weight', space=0.03)
    vis.matrix_plot(vmin=0, vmax=0.001, cbr_fontsize=9, no_colorbar=True)
    vis.plot_chromosome_bounds(linewidth=2)
    vis.plot_genes(filter_=["ABCG1"], fontsize=10, release=75)
    mkdir('./2_1_heatmap/%s' % "ABCG1")
    vis.outfig('./2_1_heatmap/%s/%s_%s.pdf' % ("ABCG1", mcool, "ABCG1"), dpi=300)

Samples = ["CRC-05-N-HiC", "CRC-05-T-HiC", "CRC-05-M-HiC"]
for sample in Samples:
    plot_heat(sample, "15,31300000,+,21,32400000,- 15,20000000 21,48100000")


