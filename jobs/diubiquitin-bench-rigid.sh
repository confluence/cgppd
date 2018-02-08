#PBS -q GPUQ
#PBS -l nodes=srvslsgpu004:ppn=20:seriesGPUk
#PBS -N 2ubq_bench_rigid
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_benchmark_rigid -t 20 -s 20 -g 2 -v 3
