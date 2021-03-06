#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=1:seriesGPUk
#PBS -N ala-1024-rep
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljrep -f config/ala1024 --gpuoffset 1 -o lj_rep/ala1024 -m 40000000 -e 8000
