#PBS -q GPUQ
#PBS -l nodes=srvslsgpu004:ppn=4:seriesGPUk
#PBS -N tetraubiquitin_lys_63
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/tetraubiquitin_lys_63 -t 4 -s 4 --gpuoffset 1
