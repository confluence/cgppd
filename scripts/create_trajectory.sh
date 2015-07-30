#!/bin/bash

cat $1/pdb/*300.0K* >> $1/trajectory.pdb
sed -i 's/END/ENDMDL/' $1/trajectory.pdb
