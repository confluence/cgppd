#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=10:seriesGPUk
#PBS -N ubq_uim_trunc
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

for configfile in config/ubq_uim_trunc/*
do
    ./cgppd -f $configfile -t 20 -s 20 -g 2 -o ubq_uim_trunc/$configfile
done