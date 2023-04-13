#!/bin/bash

# Copyright (c) 2022 Insilico Medicine

if [ -z $1 ] 
then
    echo "add number of conformations"
    exit 1
fi

if [ -z $2 ] 
then
    echo "add name of folder for save results"
    exit 1
fi

nConf=$1
folder_name=$2

make prepare_mols_info && make md

./prepare_mols_info --cfg input/config.txt --cfginp input/configInp.txt --confs $nConf --out confs.txt
time ./md --cfg input/config.txt --inp input/initAtoms.inp --cfginp input/confs.txt > log

cd "calcs/"
mkdir $folder_name
cd $folder_name
mv ../../output/calcs/* .
mv ../../input/confs.txt .
mv ../../log .
cp ../../input/configInp.txt .
cd ../..
