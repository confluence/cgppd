#PBS -q GPUQ
#PBS -l nodes=srvslsgpu001:ppn=4:seriesGPU
#PBS -N diubiquitin_lys_27
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_lys_27 -t 4 -s 4
