#PBS -q GPUQ
#PBS -l nodes=srvslsgpu003:ppn=4:seriesGPUk
#PBS -N diubiquitin

cd /home/apinska/repos/cgppd/cgppd

for configfile in config/diubiquitin_lys_48  config/diubiquitin_lys_63  config/diubiquitin_met_1
do
    ./cgppd -f configfile
done
