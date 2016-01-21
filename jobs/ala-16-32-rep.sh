#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=1:seriesGPUk
#PBS -N ala-16-32-rep
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljrep -f config/ala16 --gpuoffset 1 -o lj_rep/ala16 -m 40000000 -e 8000
./cgppd_ljrep -f config/ala32 --gpuoffset 1 -o lj_rep/ala32 -m 40000000 -e 8000
