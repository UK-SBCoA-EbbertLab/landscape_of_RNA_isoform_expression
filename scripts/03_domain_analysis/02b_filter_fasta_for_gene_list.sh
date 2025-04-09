#!/bin/bash


awk 'NR==1 {next} NR==FNR {ids[$1]; next} 
     /^>/ {header=$0; id=substr($1,2); keep=(id in ids)} 
     keep {print (header=="" ? $0 : header); header=""}' $1 $2 > $3

