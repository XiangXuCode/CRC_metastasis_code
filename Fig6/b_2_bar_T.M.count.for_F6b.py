import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import numpy as np

Samples=["02", "03", "04", "05"]
for sample in Samples:
    data = pd.read_csv('./2_3_and4_and5_all_version/2_3_CRC_%s.region.T.M.count.tsv'%sample, sep='\t')
    print(data.describe())
    region = data['region']
    Tumor = data['T']*(-1)
    Metastasis = data['M']
    plt.figure(figsize=(5, 5))
    plt.barh(region, Metastasis, color='#ED7D31', height=0.8, label='Metastasis')
    plt.barh(region, Tumor, color='#4876CA', height=0.8, label='Tumor')

    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_linewidth(1)
    ax.spines['left'].set_linewidth(False)
    plt.tick_params(width=1, labelsize=4, pad=10, length=2, left=False)
    plt.xticks(ha='center', fontproperties='arial', fontsize=8)
    plt.yticks(ha='center', fontproperties='arial', fontsize=8)
    plt.ylabel('Translocation Region', labelpad=2)
    plt.xlabel('counts', labelpad=2)
    plt.legend(('Metastasis', 'Tumor'),loc="upper center", bbox_to_anchor=(0.9, 0.9),
               prop={'size': 8}, frameon=False)
    plt.subplots_adjust(left=0.2, bottom=0.1, right=0.98, top=0.95)
    pdf = PdfPages('./2_3_and4_and5_all_version/2_5_count_version_double_bar_CRC_%s.pdf'%sample)
    pdf.savefig()
    pdf.close()
    plt.close()
