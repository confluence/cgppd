#PBS -q GPUQ
#PBS -l nodes=srvslsgpu002:ppn=4
#PBS -N 4ubq_met_1_ll
#PBS -V

cd /home/apinska/repos/cgppd --gpu-offset 2

source scripts/hex_setup.sh

./cgppd -f config/tetraubiquitin_met_1_longtail_loop -t 4 -s 4 --gpu-offset 2
