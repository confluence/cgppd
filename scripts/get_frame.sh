#!/bin/bash

FRAME=$1
CLUSTER_FILE=$2

awk "BEGIN{o=0}/frame t= $FRAME/{o=1};/ENDMDL/{o=0}{if (o){print \$0}}" $CLUSTER_FILE
