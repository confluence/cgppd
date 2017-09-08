#PBS -q GPUQ
#PBS -l nodes=srvslsgpu004:ppn=4:seriesGPUk
#PBS -N octaubiquitin_lys_11
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/octaubiquitin_lys_11 -t 4 -s 4
