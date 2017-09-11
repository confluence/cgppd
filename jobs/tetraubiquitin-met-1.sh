#PBS -q GPUQ
#PBS -l nodes=srvslsgpu004:ppn=4:seriesGPUk
#PBS -N config/tetraubiquitin_met_1
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/tetraubiquitin_met_1 -t 4 -s 4
