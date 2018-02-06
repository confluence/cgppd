#PBS -q GPUQ
#PBS -l nodes=srvslsgpu002:ppn=4:seriesGPU
#PBS -N 2ubq_bench_rigid
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

./cgppd -f config/diubiquitin_benchmark_rigid -t 4 -s 4