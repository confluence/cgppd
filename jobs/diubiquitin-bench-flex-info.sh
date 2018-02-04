#PBS -q GPUQ
#PBS -l nodes=srvslsgpu002:ppn=1:seriesGPU
#PBS -N 2ubq_bench_flex_info
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_benchmark_flex -v 3
