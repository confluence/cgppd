#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=1:seriesGPUk
#PBS -N ala-64-128-rep
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljrep -f config/ala64 --gpuoffset 1 -o lj_rep/ala64
./cgppd_ljrep -f config/ala128 --gpuoffset 1 -o lj_rep/ala128
