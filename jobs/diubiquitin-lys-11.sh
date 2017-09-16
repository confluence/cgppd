#PBS -q GPUQ
#PBS -l nodes=srvslsgpu002:ppn=4
#PBS -N 2ubq_lys_11
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_lys_11 -t 4 -s 4
#./cgppd -f config/diubiquitin_lys_11 -t 4 -s 4 --gpuoffset 2
