#PBS -q GPUQ
#PBS -l nodes=srvslsgpu002:ppn=4:seriesGPU
#PBS -N 8ubq_lys_11_ll
#PBS -V

cd /home/apinska/repos/cgppd --gpu-offset 3

source scripts/hex_setup.sh

./cgppd -f config/octaubiquitin_lys_11_longtail_loop -t 4 -s 4 --gpu-offset 3
