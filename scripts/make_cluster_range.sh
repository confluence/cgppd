#!/bin/bash

sim_dir=$1
shift

if [ ! -f "$sim_dir/trajectory.pdb" ]
then
    cat $sim_dir/pdb/*300.0K* >> "$sim_dir/trajectory.pdb"
    sed -i 's/END/ENDMDL/' "$sim_dir/trajectory.pdb"
fi

for cutoff in $@
do
    cluster_dir="$sim_dir/clusters_$cutoff"
    mkdir "$cluster_dir"
    cd "$cluster_dir"
    echo "0 0" | g_cluster -f "../trajectory.pdb" -s "../trajectory.pdb" -cutoff $cutoff -cl
    cd -
done
