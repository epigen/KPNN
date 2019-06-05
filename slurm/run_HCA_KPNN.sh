#!/bin/bash

KPNN_OUTPUTS_HCA="${KPNN_OUTPUTS}/HCA/"
mkdir $KPNN_OUTPUTS_HCA
echo $KPNN_OUTPUTS_HCA

KPNN_Replicates=1

# HCA - Bone marrow
inPathData="$KPNN_INPUTS/ica_bone_marrow_h5.h5"
inPathEdges="$KPNN_INPUTS/HCA_Edgelist.csv"
inPathYs="$KPNN_INPUTS/HCA_ClassLabels.csv"

# REAL DATA
# WITHOUT DROPOUT
label="XGXBMnoDAdj"
outPath="$KPNN_OUTPUTS_HCA/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="HCA KPNN $label" --cpus-per-task=15 --partition=longq --mem=100000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --alpha=0.05 --lambd=0.1 --disableDropoutAdjust $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

# WITH DROPOUT
label="8G8BMnoDAdj"
outPath="$KPNN_OUTPUTS_HCA/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="HCA KPNN $label" --cpus-per-task=15 --partition=longq --mem=100000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --dropOut=0.8 --dropOutGenes=0.8 --alpha=0.05 --lambd=0.1 --disableDropoutAdjust $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

label="6G6BMnoDAdj"
outPath="$KPNN_OUTPUTS_HCA/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="HCA KPNN $label" --cpus-per-task=15 --partition=longq --mem=100000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --dropOut=0.6 --dropOutGenes=0.6 --alpha=0.05 --lambd=0.1 --disableDropoutAdjust $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

label="4G4BMnoDAdj"
outPath="$KPNN_OUTPUTS_HCA/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="HCA KPNN $label" --cpus-per-task=15 --partition=longq --mem=100000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --dropOut=0.4 --dropOutGenes=0.4 --alpha=0.05 --lambd=0.1 --disableDropoutAdjust $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done





# HCA - Cord blood
inPathData="$KPNN_INPUTS/ica_cord_blood_h5.h5"
inPathEdges="$KPNN_INPUTS/HCA_Edgelist.csv"
inPathYs="$KPNN_INPUTS/HCA_ClassLabels.csv"

# REAL DATA
# WITHOUT DROPOUT
label="XGXCBnoDAdj"
outPath="$KPNN_OUTPUTS_HCA/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="HCA KPNN $label" --cpus-per-task=15 --partition=longq --mem=100000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --alpha=0.05 --lambd=0.1 --disableDropoutAdjust $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

# WITH DROPOUT
label="8G8CBnoDAdj"
outPath="$KPNN_OUTPUTS_HCA/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="HCA KPNN $label" --cpus-per-task=15 --partition=longq --mem=100000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --dropOut=0.8 --dropOutGenes=0.8 --alpha=0.05 --lambd=0.1 --disableDropoutAdjust $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

label="6G6CBnoDAdj"
outPath="$KPNN_OUTPUTS_HCA/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="HCA KPNN $label" --cpus-per-task=15 --partition=longq --mem=100000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --dropOut=0.6 --dropOutGenes=0.6 --alpha=0.05 --lambd=0.1 --disableDropoutAdjust $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done

label="4G4CBnoDAdj"
outPath="$KPNN_OUTPUTS_HCA/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="HCA KPNN $label" --cpus-per-task=15 --partition=longq --mem=100000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --dropOut=0.4 --dropOutGenes=0.4 --alpha=0.05 --lambd=0.1 --disableDropoutAdjust $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done




# CONTROL INPUTS
inPathData="$KPNN_INPUTS/ica_cord_blood_h5.h5"
inPathEdges="$KPNN_INPUTS/HCA_Edgelist.csv"
inPathYs="$KPNN_INPUTS/HCA_ClassLabels.csv"

label="Control"
outPath="$KPNN_OUTPUTS_HCA/${label}/"
mkdir $outPath
echo $outPath
for j in {1..3}; do
  sbatch --job-name="HCA KPNN $label" --cpus-per-task=10 --partition=longq --mem=70000 --time=14-00:00:00 --exclude=n[001-007] --nodes=1 \
    --wrap="for i in {1..$KPNN_Replicates}; do python $KPNN_CODEBASE/KPNN_Function.py --alpha=0.05 --lambd=0.1 --control --disableDropoutAdjust $inPathData $inPathEdges $inPathYs $outPath; done" \
    --output="$outPath/${label}_$j.log"
done
