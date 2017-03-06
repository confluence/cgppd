#PBS -q GPUQ
#PBS -l nodes=srvslsgpu004:ppn=10:seriesGPUk
#PBS -N diubiquitin-2gpus
#PBS -V

cd /home/apinska/repos/cgppd

source scripts/hex_setup.sh

for configfile in config/diubiquitin_lys_48  config/diubiquitin_lys_63  config/diubiquitin_met_1
do
    ./cgppd -f $configfile -t 20 -s 20 -g 2
done
