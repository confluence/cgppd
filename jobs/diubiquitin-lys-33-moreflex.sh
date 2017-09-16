#PBS -q GPUQ
#PBS -l nodes=srvslsgpu002:ppn=4:seriesGPU
#PBS -N 2ubq_lys_33_ll
#PBS -V

cd /home/apinska/repos/cgppd --gpu-offset 1

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_lys_33_longtail_loop -t 4 -s 4 --gpu-offset 1
