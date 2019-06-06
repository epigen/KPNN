# KPNNs
Knowledge-primed neural networks developed in the [Bock lab](http://medical-epigenomics.org) at [CeMM](http://cemm.at).

# System requirements
1. KPNNs were developed on a linux computing cluster and a Mac.

2. Training of KPNNs is performed by a python script (KPNN_Function.py). This program has been developed and tested using python 2.7.6 with the following packages: 
  - tensorflow (1.3.1)
  - tensorflow-tensorboard (0.1.8)
  - pandas (0.19.2)
  - numpy (1.13.2)
  - scipy (0.14.0)
  - tables (3.4.4)
  - psutil (5.4.6)

3. Downstream analysis is performed in R (version 3.2.3) based on the packages:
  - ggplot2 (2.2.1)
  - data.table (1.11.4)

# Installation
1. Clone this github repository
2. Install the systems requirements listed above
3. Edit the file setup.sh to define the location of the github repository (code) and where input and output data should be stored. 
4. Source setup.sh, this will create the required environmental variables (all start with "KPNN_").

# Demo to train one network
1. Download the Demo data from http://kpnn.computational-epigenetics.org/. If wget is set up on your system, you can do this using the script Download_Data.sh.
2. To train a KPNN, run the python program with four arguments: (1) Input data, (2) an edge list, (3) class labels, (4) a path to store the outputs
      ```
		  data="$KPNN_INPUTS/SIM1_Data.csv"
		  edges="$KPNN_INPUTS/SIM1_Edgelist.csv"
		  classLabels="$KPNN_INPUTS/SIM1_ClassLabels.csv"
		  outPath=$KPNN_OUTPUTS
          python $KPNN_CODEBASE/KPNN_Function.py --alpha=0.001 --lambd=0.2 $data $edges $classLabels $outPath
      ```
3. The program will produce a folder called "run_1" (or a number > 1 if another run already exists) in the location defined by $outPath, where the output of the KPNN is stored. Outputs include:
  - tf_cost.csv - tracks training progress over iterations
  - tf_NumGradMeans.csv - Node weights (from numerical gradient estimation)
  - tf_NumGradTestError.txt - Test error of the latest saved model that was used to calculate node weights
  - tf_settings.csv - The hyperparameters used to train this model
  - tf_weights.csv - Edge weights of the final network
  - tf_yHat_val.csv - Predicted class labels on test data
  - tf_yTrue_val.csv - True class labels on test data
  
# Demo to train multiple network replicates
1. To train multiple KPNNs in parallel, examples scripts to do so are provided in the folder slurm/. These scripts are based on [SLURM](slurm.schedmd.com). The script to run simulated demo networks is provided in run_SIM_KPNN.sh.
2. After training is complete, to summarize test error across trained networks, use the R script Analysis_CollectOutputs.R
3. Finally, to summarize and plot node weights in trained networks, use the R script Analysis_SIM.R (this requires running of Analysis_CollectOutputs.R first).

# Instructions on how to run KPNNs on your data
1. Example scripts to run additional datasets, such as the dataset on T cell receptor (TCR) stimulation or on predicting cell types in the Human Cell Atlas (HCA), or provided in the folder slurm/.
2. All required inputs can be downloaded using the script Download_Data.sh, except for single-cell expression data from the Human Cell Atlas, which should be downloaded under https://preview.data.humancellatlas.org and should then be stored under $KPNN_INPUTS (defined in setup.sh)
3. After training, use Analysis_CollectOutputs.R, and then Analysis_TCR.R or Analysis_HCA.R to summarize the results across network replicates.
4. To adjust these analyses to your dataset, you mainly need to adjust the inputs to provide your data, class labels, and edgelist.
