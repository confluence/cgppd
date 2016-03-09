#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=6:seriesGPUk
#PBS -N ala-128-ljrep-replica-exchange
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljrep -f config/ala128 -t 20 -s 20 -o lj_rep/ala128 -m 40000000 -e 8000 -n 240 -x 420
