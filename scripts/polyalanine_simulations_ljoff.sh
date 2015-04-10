#!/usr/bin/bash

make clean
make LJ=off STREAMS=yes

for n in 4 8 16 32 64 128 256 512 1024 2048
do
    ./cgppd_ljoff -f config/ala$n
done

# TODO: add processing (pdbstats?)
