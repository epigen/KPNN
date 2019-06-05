#!/bin/bash

# TCR
KPNN_OUTPUTS_TCR="${KPNN_OUTPUTS}/TCR/"
mkdir $KPNN_OUTPUTS_TCR
echo $KPNN_OUTPUTS_TCR

KPNN_Replicates=2

# HCA - Bone marrow
inPathData="$KPNN_INPUTS/TCR_Data.h5"
inPathEdges="$KPNN_INPUTS/TCR_Edgelist.csv"
inPathYs="$KPNN_INPUTS/TCR_ClassLabels.csv"

# REAL DATA
# WITHOUT DROPOUT
label="XGX"
outPath="$KPNN_OUTPUTS_TCR/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="TCR $label" --cpus-per-task=10 --partition=longq --mem=60000 --time=14-00:00:00 --exclude=i[001-020],n[001-005] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --alpha=0.01 --lambd=0.1 --maxBreakCount=20 --minImprovement=0.8 $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

# WITH DROPOUT
ARR=( 4 5 6 7 8 9 )
for d in ${ARR[@]}; do
  label="${d}G${d}"
  outPath="$KPNN_OUTPUTS_TCR/${label}/"
  mkdir $outPath
  echo $outPath
  for j in {1..3}; do
    sbatch --job-name="TCR $label" --cpus-per-task=10 --partition=longq --mem=60000 --time=14-00:00:00 --exclude=i[001-020],n[001-005] --nodes=1 \
      --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --dropOut=0.${d} --dropOutGenes=0.${d} --alpha=0.01 --lambd=0.1 --maxBreakCount=20 --minImprovement=0.8 $inPathData $inPathEdges $inPathYs $outPath; done" \
      --output="$outPath/${label}_$j.log"
  done
done


# CONTROL INPUTS
# WITHOUT DROPOUT
label="XGXControl"
outPath="$KPNN_OUTPUTS_TCR/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="TCR $label" --cpus-per-task=10 --partition=longq --mem=60000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --alpha=0.01 --lambd=0.1 --maxBreakCount=20 --minImprovement=0.8 --control $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

# WITH DROPOUT
ARR=( 4 5 6 7 8 9 )
for d in ${ARR[@]}; do
  label="${d}G${d}Control"
  outPath="$KPNN_OUTPUTS_TCR/${label}/"
  mkdir $outPath
  echo $outPath
  for j in {1..3}; do
    sbatch --job-name="TCR $label" --cpus-per-task=10 --partition=longq --mem=60000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
      --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --dropOut=0.${d} --dropOutGenes=0.${d} --alpha=0.01 --lambd=0.1 --maxBreakCount=20 --minImprovement=0.8 --control $inPathData $inPathEdges $inPathYs $outPath; done" \
      --output="$outPath/${label}_$j.log"
  done
done

