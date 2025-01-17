#!/bin/bash
#SBATCH -J apacp
#SBATCH -p gpu_4l
#SBATCH -N 1 
#SBATCH -o apa_cp1loops_%j.out
#SBATCH -e app_cp1loops_%j.err
#SBATCH --no-requeue
#SBATCH -A lch3000_g1
#SBATCH --qos=lch3000gpu
#SBATCH -c 1
hic_path="/lustre3/lch3000_pkuhpc/ganjb/Shougang/loop"
loop_path="/lustre3/lch3000_pkuhpc/ganjb/Shougang/loop"

time java -Xmx8g -jar lustre3/lch3000_pkuhpc/ganjb/software/juicer_tools_1.11.09_jcuda.0.8.jar apa  \-r 10000 -n 10 -w 5 -u  \
${hic_path}/all_CN.validPairs.hic \${loop_path}/CN/merged_loops.bedpe \${loop_path}/CN/APA/CN-allloops-n10-w5

###this script is used to calculate APA score for selected loops on cluster (GPU is required)
