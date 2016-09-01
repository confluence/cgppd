#!/bin/bash

if [ -z "$2" ]
then
    temp=300.0
else
    temp=$2
fi

if test -n "$(find $1/pdb -maxdepth 1 -name *${temp}K* -print -quit)"
then
    rm -f $1/trajectory.pdb
    cat $1/pdb/*${temp}K* >> $1/trajectory.pdb
    sed -i 's/END/ENDMDL/' $1/trajectory.pdb
else
    echo "No samples found at temperature ${temp}K."
fi





