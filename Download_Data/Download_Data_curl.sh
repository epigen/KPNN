#!/bin/bash

web_input=$(curl -Ls -o /dev/null -w %{url_effective} http://kpnn.computational-epigenetics.org/)
echo $web_input

# TCR DATA
ARR=( TCR_Edgelist.csv TCR_ClassLabels.csv TCR_Data.h5 HCA_ClassLabels.csv HCA_Edgelist.csv SIM1_ClassLabels.csv SIM1_Data.csv SIM1_Edgelist.csv SIM2_ClassLabels.csv SIM2_Data.csv SIM2_Edgelist.csv)
for file in ${ARR[@]}; do
    curl $web_input/$file -o $KPNN_INPUTS/$file
done