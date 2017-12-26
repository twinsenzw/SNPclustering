#!/bin/bash

cd /home/zhouw/sepi/SNP_mapping/strains

#----1. align strain genomes, identify core blocks, collect SNP information----

source /home/zhouw/mugsy/mugsyenv.sh
/home/zhouw/mugsy/mugsy --directory /home/zhouw/sepi/SNP_mapping/mugsy --prefix alignment *.fa




