#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=1:seriesGPUk
#PBS -N ala-4-8-off
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljoff -f config/ala4 -o lj_off/ala4
./cgppd_ljoff -f config/ala8 -o lj_off/ala8
