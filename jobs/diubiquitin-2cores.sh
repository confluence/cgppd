#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=2:seriesGPUk
#PBS -N diubiquitin
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

for configfile in config/diubiquitin_lys_48  config/diubiquitin_lys_63  config/diubiquitin_met_1
do
    ./cgppd -f $configfile -t 2 -s 2 --gpuoffset 1
done