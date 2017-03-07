#PBS -q GPUQ
#PBS -l nodes=srvslsgpu004:ppn=4:seriesGPUk
#PBS -N diubiquitin_lys_11
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_lys_11 -t 4 -s 4
#./cgppd -f config/diubiquitin_lys_11 -t 4 -s 4 --gpuoffset 2
