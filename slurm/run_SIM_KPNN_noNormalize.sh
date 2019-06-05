#!/bin/bash


#
# SIM1
#
KPNN_OUTPUTS_SIM="${KPNN_OUTPUTS}/SIM1N/"
mkdir $KPNN_OUTPUTS_SIM
echo $KPNN_OUTPUTS_SIM

# HCA - Bone marrow
inPathData="$KPNN_INPUTS/SIM1_Data.csv"
inPathEdges="$KPNN_INPUTS/SIM1_Edgelist.csv"
inPathYs="$KPNN_INPUTS/SIM1_ClassLabels.csv"

label="SIM1"
outPath="$KPNN_OUTPUTS_SIM/${label}/"
mkdir $outPath
for j in {1..3}; do
sbatch --job-name="SIM KPNN $label" --cpus-per-task=10 --partition=longq --mem=60000 --time=7-00:00:00 --exclude=i[001-020] --nodes=1 \
    --wrap="for i in {1..30}; do python $KPNN_CODEBASE/KPNN_Function.py --alpha=0.001 --lambd=0.2 --normalizedMatrix --normalizedMatrix01 $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done


#
# SIM2
#
KPNN_OUTPUTS_SIM="${KPNN_OUTPUTS}/SIM2N/"
mkdir $KPNN_OUTPUTS_SIM
echo $KPNN_OUTPUTS_SIM

# HCA - Bone marrow
inPathData="$KPNN_INPUTS/SIM2_Data.csv"
inPathEdges="$KPNN_INPUTS/SIM2_Edgelist.csv"
inPathYs="$KPNN_INPUTS/SIM2_ClassLabels.csv"

label="SIM2NoDropout"
outPath="$KPNN_OUTPUTS_SIM/${label}/"
mkdir $outPath
for j in {1..3}; do
sbatch --job-name="SIM KPNN $label" --cpus-per-task=10 --partition=longq --mem=60000 --time=7-00:00:00 --exclude=i[001-020] --nodes=1 \
    --wrap="for i in {1..30}; do python $KPNN_CODEBASE/KPNN_Function.py --alpha=0.05 --lambd=0.2 --normalizedMatrix --normalizedMatrix01 $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

label="SIM2Dropout"
outPath="$KPNN_OUTPUTS_SIM/${label}/"
mkdir $outPath
for j in {1..3}; do
sbatch --job-name="SIM KPNN $label" --cpus-per-task=10 --partition=longq --mem=60000 --time=7-00:00:00 --exclude=i[001-020] --nodes=1 \
    --wrap="for i in {1..30}; do python $KPNN_CODEBASE/KPNN_Function.py --alpha=0.05 --lambd=0.2 --normalizedMatrix --normalizedMatrix01 --dropOut=0.7 $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

