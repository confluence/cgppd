#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=4:seriesGPUk
#PBS -N config/diubiquitin_met_1
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_met_1 -t 4 -s 4 --gpuoffset 2
