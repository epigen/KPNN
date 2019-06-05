#!/bin/bash

web_input=$(curl -Ls -o /dev/null -w %{url_effective} http://kpnn.computational-epigenetics.org/)
echo $web_input

# TCR DATA
ARR=( TCR_Edgelist.csv TCR_ClassLabels.csv TCR_Data.h5 HCA_ClassLabels.csv HCA_Edgelist.csv )
for file in ${ARR[@]}; do
    wget $web_input/$file -O $KPNN_INPUTS/$file
done