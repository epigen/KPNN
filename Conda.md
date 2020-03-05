# Conda
```
	conda create -n kpnn python=3.7 tensorflow=1.13 pandas=0.24 numpy=1.16 scipy=1.2 psutil=5.6
	conda activate kpnn
	conda install -c conda-forge pytables=3.5
	conda install -c r r
	source Demo_setup.sh
	Rscript Requirements_R.R
	Rscript Test_Data.R $KPNN_INPUTS
	python KPNN_Function.py
```