#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=10:seriesGPUk
#PBS -N 2ubq_bench_half
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_benchmark_half -t 10 -s 10 -v 3 -g 4
