#!/bin/bash

sed -ri 's/(ATOM .{49}).{6}(.*)/\1  0.00\2/' $1

for res in `seq $2 $3`
do
    sed -ri "s/(ATOM .{17} {0,3}$res .{27}).{6}(.*)/\1  1.00\2/" $1
done
