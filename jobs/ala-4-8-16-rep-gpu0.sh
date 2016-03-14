#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=6:seriesGPUk
#PBS -N ala-4-8-16-ljrep-replica-exchange
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd_ljrep -f config/ala4 -r 20 -t 20 -s 20 -o lj_rep/ala4 -m 40000000 -e 8000 -n 240 -x 420
./cgppd_ljrep -f config/ala8 -r 20 -t 20 -s 20 -o lj_rep/ala8 -m 40000000 -e 8000 -n 240 -x 420
./cgppd_ljrep -f config/ala16 -r 20 -t 20 -s 20 -o lj_rep/ala16 -m 40000000 -e 8000 -n 240 -x 420
