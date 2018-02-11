#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=12:seriesGPUk
#PBS -N 2ubq_bench_rigid
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_benchmark_rigid -t 12 -s 12 -v 3 -g 4
