#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=1:seriesGPUk
#PBS -N ala-1024-off
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljoff -f config/ala1024
