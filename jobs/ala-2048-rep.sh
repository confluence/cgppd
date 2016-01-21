#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=1:seriesGPUk
#PBS -N ala-2048-rep
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljrep -f config/ala2048 --gpuoffset 1 -o lj_rep/ala2048 -m 40000000 -e 8000
