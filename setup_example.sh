#!/bin/bash

# Point to the base of the github repository
export KPNN_CODEBASE=$PWD

# Where input data are saved (expression data, edgelists, class labels)
export KPNN_INPUTS=$PWD/Inputs/         
mkdir -p $KPNN_INPUTS

# Where outputs are saved (trained KPNNs)
export KPNN_OUTPUTS=$PWD/Outputs/       
mkdir -p $KPNN_OUTPUTS

# Results of demo runs
export KPNN_DEMOS=$PWD/Demos/           
mkdir -p $KPNN_DEMOS

# Results of Analyses of KPNN outputs
export KPNN_ANALYSES=$PWD/Analysis/     
mkdir -p $KPNN_ANALYSES