# KPNNs
Knowledge-primed neural networks developed in the [Bock lab](http://medical-epigenomics.org) at [CeMM](http://cemm.at). KPNNs are neural networks that learn to predict class values from input data given a specific network structure. After training, they enable extraction of extraction of node weights using numerical gradient estimation. Node weights represent the importance of individual nodes for predictions.


# System requirements
1. KPNNs were developed on linux and on Mac.
2. Training of KPNNs is performed by a python script (KPNN_Function.py). This program has been developed and tested using python (tested on versions 2.7.6 and 2.7.13) 
3. Downstream analysis is performed in R (tested on versions 3.2.3 or 3.5.1)

# Installation
1. Make sure you use the right versions of python and R:
      ```
	  python --version
	  R --version
      ```
2. Create a directory to contain all analyses:
      ```
	  mkdir KPNN/
	  cd KPNN/
      ```
3. The recommended approach is to use a virtual environment:
      ```
	  # pip install virtualenv # if not installed yet
	  mkdir virtenv
	  virtualenv virtenv --no-site-packages
	  source virtenv/bin/activate
      ```
4. Installation instructions:
	  ```
	  # Clone this repository
	  git clone https://github.com/epigen/KPNN.git
	  cd KPNN/
	  
	  # Sets up environmental variables
	  source setup_example.sh

	  # Install requirements (python and R)
	  pip install -r Requirements_python.txt
	  Rscript Requirements_R.R

	  # Generate test data and run test
	  Rscript Generate_test_data.R
	  python KPNN_Function.py
      ```
5. The last command should finish without errors and you should see the message "KPNN TRAINING COMPLETED SUCCESSFULLY". Then installation is complete.

# Demo to train one network
1. Download the Demo data from http://kpnn.computational-epigenetics.org/. If wget or curl is set up on your system, you can do this using the scripts under Download_Data/.
2. To train a KPNN, run the python program with four arguments: (1) Input data, (2) an edge list, (3) class labels, (4) a path to store the outputs, for example:
      ```
		  data="$KPNN_INPUTS/SIM1_Data.csv"
		  edges="$KPNN_INPUTS/SIM1_Edgelist.csv"
		  classLabels="$KPNN_INPUTS/SIM1_ClassLabels.csv"
		  outPath=$KPNN_DEMOS
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
4. The program should run for approximately 30 minutes or less
  
# Demo to train multiple network replicates
1. To train multiple KPNNs in parallel, examples scripts to do so are provided in the folder Slurm/. These scripts are based on [SLURM](slurm.schedmd.com). The script to run simulated demo networks is provided in run_SIM_KPNN.sh.
2. After training is complete, to summarize test error across trained networks, use the R script Analysis_Scripts/Analysis_CollectOutputs.R
3. Finally, to summarize and plot node weights in trained networks, use the R script Analysis_Scripts/Analysis_SIM.R (this requires running of Analysis_CollectOutputs.R first).

# Instructions on how to run KPNNs on your data
1. Example scripts to run additional datasets, such as the dataset on T cell receptor (TCR) stimulation or on predicting cell types in the Human Cell Atlas (HCA), or provided in the folder Slurm/.
2. All required inputs can be downloaded using the scripts in Download_Data/, except for single-cell expression data from the Human Cell Atlas, which should be downloaded under https://preview.data.humancellatlas.org and should then be stored under $KPNN_INPUTS (defined in setup.sh)
3. After training, use scripts in Analysis_Scripts such as Analysis_CollectOutputs.R, and then Analysis_TCR.R or Analysis_HCA.R to summarize the results across network replicates.
4. To adjust these analyses to your dataset, you need to adjust the inputs to provide your data, class labels, and edgelist.
