#!/bin/bash

web_input=$(curl -Ls -o /dev/null -w %{url_effective} http://kpnn.computational-epigenetics.org/)
echo $web_input

# TCR DATA
ARR=(SIM1_ClassLabels.csv SIM1_Data.csv SIM1_Edgelist.csv)
for file in ${ARR[@]}; do
    curl $web_input/$file -o $KPNN_INPUTS/$file
done