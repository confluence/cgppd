#!/bin/bash

for l in `find $1 -name log`
do
    d=`dirname "$l"`
    tar -C "$d" -czf "$l.tgz" --remove-files log
done
