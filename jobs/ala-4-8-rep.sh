#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=1:seriesGPUk
#PBS -N ala-4-8-rep
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljrep -f config/ala4 --gpuoffset 1 -o lj_rep/ala4 -m 20000000 -e 4000
./cgppd_ljrep -f config/ala8 --gpuoffset 1 -o lj_rep/ala8 -m 20000000 -e 4000
