#PBS -q GPUQ
#PBS -l nodes=srvslsgpu002:ppn=4
#PBS -N 4ubq_lys_48
#PBS -V

cd /home/apinska/repos/cgppd --gpu-offset 2

source scripts/hex_setup.sh

./cgppd -f config/tetraubiquitin_lys_48 -t 4 -s 4 --gpu-offset 2
