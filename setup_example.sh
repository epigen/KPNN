#!/bin/bash

# Point to where you want input and output data to go
export KPNN_ROOT=$PWD
# Point to the base of the github repository
export KPNN_CODEBASE=$PWD

# These will be populated automatically but can also be adjusted
export KPNN_OUTPUTS=$KPNN_ROOT/Outputs/
export KPNN_INPUTS=$KPNN_ROOT/Inputs/
export KPNN_DEMOS=$KPNN_ROOT/Demos/

mkdir -p $KPNN_ROOT
mkdir -p $KPNN_OUTPUTS
mkdir -p $KPNN_INPUTS
mkdir -p $KPNN_DEMOS