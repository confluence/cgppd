#!/bin/bash

sim_dir=$1
shift

trajectory_file=$sim_dir/trajectory.pdb

if [ ! -f "$trajectory_file" ]; then
    cat $sim_dir/pdb/*300.0K* >> "$trajectory_file"
    sed -i 's/END/ENDMDL/' "$trajectory_file"
fi

for cutoff in $@
do
    cluster_dir="$sim_dir/clusters_$cutoff"
    mkdir "$cluster_dir"
    cd "$cluster_dir"
    g_cluster -f "$trajectory_file" -s "$trajectory_file" -cutoff $cutoff -cl
    cd -
done
