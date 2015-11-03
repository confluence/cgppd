#!/bin/bash

make clean
make LJ=repulsive

# This prints compilation options
./cgppd_ljrep

for n in 4 8 16 32 64 128 256 512 1024 2048
do
    ./cgppd_ljrep -f config/ala$n
done

