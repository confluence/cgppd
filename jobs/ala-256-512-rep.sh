#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=1:seriesGPUk
#PBS -N ala-256-512-rep
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljrep -f config/ala256 --gpuoffset 1 -o lj_rep/ala256 -m 20000000 -e 4000
./cgppd_ljrep -f config/ala512 --gpuoffset 1 -o lj_rep/ala512 -m 20000000 -e 4000
