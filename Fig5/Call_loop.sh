#!/bin/bash
#SBATCH -J cancer_loop
#SBATCH -p gpu_4l
#SBATCH -N 1
#SBATCH -o loop_%j.out
#SBATCH -e loop_%j.err
#SBATCH --no-requeue
#SBATCH -A lch3000_g1
#SBATCH --qos=lch3000gpu
#SBATCH --mincpus=14
#source /apps/source/gcc-4.7.4fortran.sh
#source /apps/source/cuda-8.0.61.sh
#sof="time java -Xmx10g -jar /lustre1/lch3000_pkuhpc/liuyt/software/bin/juicer_tools.1.7.6_jcuda.0.8.jar hiccups -m 1024"
sof="time java -jar /lustre3/lch3000_pkuhpc/ganjb/software/juicer_tools_1.11.09_jcuda.0.8.jar hiccups -m 1024"
par="-r 10000 \
-f 0.15 \
-p 2 \
-i 5 \
-t 0.02,1.5,1.75,2 \
-d 20000"

time ${sof} ${par} /lustre3/lch3000_pkuhpc/ganjb/Shougang/loop/all_CN.validPairs.hic /lustre3/lch3000_pkuhpc/ganjb/Shougang/loop/CN

time ${sof} ${par} /lustre3/lch3000_pkuhpc/ganjb/Shougang/loop/all_CT.validPairs.hic /lustre3/lch3000_pkuhpc/ganjb/Shougang/loop/CT

time ${sof} ${par} /lustre3/lch3000_pkuhpc/ganjb/Shougang/loop/all_LT.validPairs.hic /lustre3/lch3000_pkuhpc/ganjb/Shougang/loop/LT

###this script is used to call loop on cluster (GPU is required)
