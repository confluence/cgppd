#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=1:seriesGPUk
#PBS -N ala-256-512-off
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljoff -f config/ala256 -o lj_off/ala256
./cgppd_ljoff -f config/ala512 -o lj_off/ala512
