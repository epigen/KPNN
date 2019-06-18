import tensorflow as tf
import pandas as pd
import numpy as np
import re
import sys
import time
import io
import os
import argparse
import random as rd
import copy
import collections
import scipy.sparse as sp_sparse
import tables
import psutil
import gc



######################
## PARSE ARGUMENTS ###
######################
parser = argparse.ArgumentParser(description='Neural network.')

# Inputs
parser.add_argument('inPathData', type=str, help='path to data file')
parser.add_argument('inPathEdges', type=str, help='path to edges file')
parser.add_argument('inPathYs', type=str, help='path to output (ys) file')
parser.add_argument('outPath', type=str, help='path to output folder (must exist)')

# Data / Normalization
parser.add_argument('--genome', type=str, help='genome in hdf5 file, is guessed by default', default="")
parser.add_argument('--normalizedMatrix', action='store_false', help="do TPM, log of each cell")
parser.add_argument('--normalizedMatrix01', action='store_false', help="do 0-1 normalization of each gene")
parser.add_argument('--maxCells', type=float, help="How many cells can be run at once", default=100000)

# Training parameters
parser.add_argument('--iterations', type=int, help='number of iterations before automatic quit', default=100000)
parser.add_argument('--alpha', type=float, help='learning rate alpha', default=0.1)
parser.add_argument('--momentum', type=float, help='If set uses Momentum Optimzer in Tensorflow', default=0)
parser.add_argument('--lambd', type=float, help='regularization lambda for l2 loss', default=0.0)
parser.add_argument('--testSet', type=float, help='fraction of samples for test set', default=0.2)
parser.add_argument('--minibatch', type=int, help="Size of minibatch", default=1000)

# Early stopping
parser.add_argument('--disableInterrupt', action='store_true', help="when toggled removes early stopping")
parser.add_argument('--maxBreakCount', type=int, help='maximum number of non succesful runs before we stop training', default=10)
parser.add_argument('--minImprovement', type=float, help='fraction improvement over best run to save model', default=0.7)

# Dropout
parser.add_argument('--dropOut', type=float, help="Dropout Keep Probability Nodes", default=1.0)
parser.add_argument('--dropOutGenes', type=float, help="Dropout Keep Probability Genes", default=1.0)
parser.add_argument('--disableDropoutAdjust', action='store_true', help="when toggled dropout is not adjusted for child nodes of parent nodes with few childen")

# Control / shuffle /set seed /
parser.add_argument('--shuffleGenes', action='store_true', help="Shuffle gene expression")
parser.add_argument('--control', action='store_true', help="Run on control data (all genes are predictive)")
parser.add_argument('--randomSeed', type=int, help='Seed to set for numpy and tensorflow (this only works in python2 for some reason)')

# Tracking
parser.add_argument('--tfWrite', action='store_true', help="Write train and test summaries")

# Numerical gradient estimation of node importance
parser.add_argument('--numGradEpsilon', type=float, help="Epsilon of numerical gradient approximation, done on all nodes in the network", default=0.001)
parser.add_argument('--disableNumGrad', action='store_true', help="Do not perform numerical gradient estimation of node importance (saves quite some time)")

# Computing limits
parser.add_argument('--threads', type=int, help="Parallelization", default=10)



##############################################
## DEFAULT ARGUMENTS FOR TEST RUNS ###########
##############################################
if len(sys.argv) < 4: # means we are in python shell or the script is run without arguments
    args = parser.parse_args([
        os.environ['KPNN_INPUTS'] + "/TEST_Data.csv",
        os.environ['KPNN_INPUTS'] + "/TEST_Edgelist.csv",
        os.environ['KPNN_INPUTS'] + "/TEST_ClassLabels.csv",
        os.environ['TMPDIR']
    ])
    args.lambd = 0.01
    args.iterations = 5
    args.threads = 1
    args.control = False
    args.randomSeed=1
else:           # script is being called from outside with proper arguments
    args = parser.parse_args()

print(args)



##############################################
## RANDOM SEED ###########
##############################################
if(args.randomSeed is not None):
    print("----------------\n...Setting random seed")
    np.random.seed(args.randomSeed)
    tf.random.set_random_seed(args.randomSeed)
    
    #Test random seeds
    sess2 = tf.Session()
    weights = tf.Variable(tf.random_normal([1,1], dtype=tf.float64, name="Random_weights"), name="Weights",dtype=tf.float64)
    init = tf.global_variables_initializer()
    sess2.run(init)
    print("Tensorflow random number:" + str(sess2.run(weights)))
    print("Numpy choice:")
    print(np.random.choice(list(range(1000))))
    print("Numpy Binomial")
    print(np.random.binomial(list(range(1000)), 0.4)[200])
    print("Numpy Shuffle")
    x = list(range(1000))
    np.random.shuffle(x)
    print(x[5])


##############################
## SETUP ARGUMENTS ###########
##############################
print("----------------\n...Processing arguments")
# ADAM
doADAM = True
if args.momentum > 0: doADAM = False

# Dropout on nodes
doDropout = args.dropOut < 1
if doDropout: 
    print("Using dropout on Nodes " + str(args.dropOut))

# Dropout on genes
doDropoutGenes = args.dropOutGenes < 1
if doDropoutGenes: 
    print("Using dropout on Genes " + str(args.dropOutGenes))

# Figure out where to put output
if not os.path.exists(args.outPath):
    os.mkdir(args.outPath)

i = 1
while os.path.exists(os.path.join(args.outPath, "run_" + str(i))):
    i += 1

outPath = os.path.join(args.outPath, "run_" + str(i))
os.mkdir(outPath)



##############################
## TRACK MEMORY FUNCTION #####
##############################
open(os.path.join(outPath, "Memory.tsv"),"w").write("Step\tMemory_MB\tTime\n")
def logMem(text):
    open(os.path.join(outPath, "Memory.tsv"),"a").write(text + "\t" + str(psutil.Process(os.getpid()).memory_info().rss*1.0/10**6) + "\t" + time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + "\n")

logMem("Start")



##############################
## LOAD ALL DATA #############
##############################
print("----------------\n...Loading Data")

# SET UP FUNCTIONS TO LOAD DATA --------------------------------------------------------------------------------------------------------------------
GeneBCMatrix = collections.namedtuple('GeneBCMatrix', ['gene_ids', 'gene_names', 'barcodes', 'matrix'])
def get_matrix_from_h5(filename, genome):
    with tables.open_file(filename, 'r') as f:
        try:
            group = f.get_node(f.root, genome)
        except tables.NoSuchNodeError:
            print("That genome does not exist in this file.")
            return None
        gene_ids = getattr(group, 'genes').read()
        gene_names = getattr(group, 'gene_names').read()
        barcodes = getattr(group, 'barcodes').read()
        data = getattr(group, 'data').read()
        indices = getattr(group, 'indices').read()
        indptr = getattr(group, 'indptr').read()
        shape = getattr(group, 'shape').read()
        matrix = sp_sparse.csc_matrix((data, indices, indptr), shape=shape)
        return GeneBCMatrix(gene_ids, gene_names, barcodes, matrix)

def indexInList(l2, l1):
    assert len(set(l1)) == len(l1), "Provided list contains duplicated elements"
    ref={}
    res=[]
    for idx, x in enumerate(l1):
        ref[x]=idx
    
    for idx, x in enumerate(l2):
        res.append(ref[x])
    
    return(res)

# LOAD INPUT DATA --------------------------------------------------------------------------------------------------------------------
if re.match(".+csv$", args.inPathData) is not None:
    df=pd.read_csv(args.inPathData, sep=",", index_col=0)
    genesList_x = df.index.tolist()
    assert all([genesList_x[idx] == df.index[idx] for idx in range(len(genesList_x))]) # Geneslist corresponds to column names
    assert len(genesList_x) == len(set(genesList_x)) # they are unique
    barcodes_x = df.columns.tolist()
    fullData = sp_sparse.csc_matrix(df.as_matrix().astype("float64"))
elif re.match(".+h5$", args.inPathData) is not None:
    if args.genome == "":
        h5f = tables.open_file(args.inPathData, 'r')
        args.genome = h5f.list_nodes(h5f.root)[0]._v_name
        h5f.close()
    gene_bc_matrix = get_matrix_from_h5(args.inPathData, args.genome)
    barcodes_x = gene_bc_matrix.barcodes.tolist()
    genesList_x = gene_bc_matrix.gene_names.tolist()
    fullData = gene_bc_matrix.matrix.astype("float64")
else:
    assert False, "No adequate input data provided"

genesList_x = [x + "_gene" for x in genesList_x]
logMem("Data loaded")

# LOAD Ys and match barcodes to input data --------------------------------------------------------------------------------------------------------------------
file_y=pd.read_csv(args.inPathYs, sep=",")
barcodes_y = file_y["barcode"].tolist()
# print(barcodes_y[0:5])
# print(set(barcodes_y[0:5]))
# print(barcodes_x[0:5])
# print(set(barcodes_x[0:5]))
# match barcodes:
barcodes = sorted(set(barcodes_y) & set(barcodes_x))
# print(barcodes[0:5])
assert len(barcodes) == len(set(barcodes)), "BARCODES ARE NOT UNIQUE!"
print("Number of Barcodes found in Y but not X: " + str(len(set(barcodes_y) - set(barcodes_x))))
print("Number of Barcodes used: " + str(len(barcodes)))
file_y=file_y.loc[indexInList(barcodes, barcodes_y)]
open(os.path.join(outPath, "barcodes.txt"),"w").write("\n".join(barcodes))

# LOAD EDGE LIST and match gene list with input data --------------------------------------------------------------------------------------------------------------------
edgelistFile=pd.read_csv(args.inPathEdges, sep=",")
genesList_edges = list(set(edgelistFile['child'].tolist()) - set(edgelistFile['parent'].tolist()))
genesList = sorted(set(genesList_x) & set(genesList_edges))
print("Number of Genes used: " + str(len(genesList)))
assert len(genesList) == len(set(genesList)), "genesList non unique!"
genesListOrig = copy.copy(genesList)
edgelistFile = edgelistFile[edgelistFile['child'].isin(genesList) | edgelistFile['child'].isin(edgelistFile['parent'])]
open(os.path.join(outPath, "genesList.txt"),"w").write("\n".join(genesList))

# NETWORK OUTPUTS matched between Ys and edgelist --------------------------------------------------------------------------------------------------------------------
outputs_y = sorted(set(file_y.columns.tolist()) - set(["barcode"]))
outputs_edges = sorted(set(edgelistFile['parent'].tolist()) - set(edgelistFile['child'].tolist()))
outputs = sorted(set(outputs_y) & set(outputs_edges))
print("Number of Outputs used: " + str(len(outputs))  + " --> " + ",".join(outputs))
print("\tClass labels: " + ",".join(outputs_y))
print("\tNetwork outputs: " + ",".join(outputs_edges))
assert len(outputs) > 0, "No outputs fitting between y table and edgelist"
assert len(outputs) == len(set(outputs)), "outputs non unique!"
open(os.path.join(outPath, "outputs.txt"),"w").write("\n".join(outputs))

# MAKE SURE EDGES MAKE UP DAG --------------------------------------------------------------------------------------------------------------------
outputs_edges_missed = list(set(outputs_edges) - set(outputs))
input_edges_missed = list(set(set(edgelistFile['child'].tolist()) - set(edgelistFile['parent'].tolist())) - set(genesList))
while len(outputs_edges_missed) > 0 or len(input_edges_missed) > 0:
    edgelistFile = edgelistFile[~edgelistFile['parent'].isin(outputs_edges_missed)]
    edgelistFile = edgelistFile[~edgelistFile['child'].isin(input_edges_missed)]
    # These are the parents that are not children and also not in outputs (leaves with only out edges - which can only be outputs)
    outputs_edges_missed = list(set(set(edgelistFile['parent'].tolist()) - set(edgelistFile['child'].tolist())) - set(outputs))
    # These are the children that are not parents and also not in the gene list (leaves with only in edges - which can only be genes)
    input_edges_missed = list(set(set(edgelistFile['child'].tolist()) - set(edgelistFile['parent'].tolist())) - set(genesList))

for x in outputs:
    assert x in edgelistFile['parent'].tolist(), x + " missing from edgelist"

edgelistFile.to_csv(os.path.join(outPath,"edges.tsv"), index=False)

# TRANFORM Ys to fit to our data --------------------------------------------------------------------------------------------------------------------
fullY = np.transpose(file_y[outputs].as_matrix()).astype("float64")

# GET RELEVANT INPUT DATA FROM FULL DATASET --------------------------------------------------------------------------------------------------------------------
fullData = fullData[:,indexInList(barcodes, barcodes_x)]

# GENERATE CONTROL DATA --------------------------------------------------------------------------------------------------------------------
if args.control:
    # edgelist
    edgelistFile = edgelistFile.loc[~edgelistFile['parent'].isin(outputs[1:])]
    edgelistFile.parent[edgelistFile['parent'].isin(outputs[0:1])] = "output"
    outputs = ["output"]
    
    # Set up fullY
    size_ds = int(min(fullY.shape[1], args.minibatch*2)/2) * 2 # this must be devidable by two
    fullY = fullY[0:1,:size_ds]
    fullY[:,:size_ds//2] = 1
    fullY[:,size_ds//2:] = 0
    
    # set up Barcodes
    barcodes = barcodes[:size_ds]
    barcodes_x = barcodes
    
    # Sum up expression to draw from
    control_grp1 = np.asarray(fullData.sum(1))[:,0].astype("int64")
    control_grp2 = control_grp1 * np.random.choice([0.5, 2], control_grp1.shape[0]).astype("int64") # add two fold up down regulation everywhere
    
    # Generate full data
    control_grp1_sum = control_grp1.sum()
    control_grp2_sum = control_grp2.sum()
    fullData_colSums = fullData.sum(0)
    fullData_control = np.zeros((control_grp1.shape[0],size_ds), dtype=int)
    assert fullData_control.shape[1] <= fullData.shape[1]
    for i2 in range(size_ds):
        if int(fullY[0,i2]) == 1:
            fullData_control[:,i2] = np.random.binomial(control_grp1, fullData_colSums[0,i2]/control_grp1_sum)
        else:
            fullData_control[:,i2] = np.random.binomial(control_grp2, fullData_colSums[0,i2]/control_grp2_sum)
    
    fullData = sp_sparse.csc_matrix(fullData_control.astype("float64"))

# Normalize log TPM (by columms) --------------------------------------------------------------------------------------------------------------------
if args.normalizedMatrix:
    for i in range(fullData.shape[1]):
        x=fullData.data[fullData.indptr[i]:fullData.indptr[i+1]]
        fullData.data[fullData.indptr[i]:fullData.indptr[i+1]]=x/x.sum()*1e6
    fullData = fullData.log1p()

# Extract relevant genes --------------------------------------------------------------------------------------------------------------------
fullData = fullData[indexInList(genesList, genesList_x),:]

# Define train and test data  --------------------------------------------------------------------------------------------------------------------
test_idx = []
val_idx = []
train_idx = []
# If we have defined train/test/validation set in the y_file then use those instead
if "Set" in list(file_y): # if some elements are defined
    assert file_y["barcode"].tolist() == barcodes
    if set(file_y['Set'].tolist()) == set(['train', 'test', 'val']): # if all sets are defined
        print("Setting train - val - test sets from y-file")
        test_idx = indexInList(file_y[file_y['Set']=="test"]['barcode'].tolist(), barcodes)
        val_idx = indexInList(file_y[file_y['Set']=="val"]['barcode'].tolist(), barcodes)
        train_idx = indexInList(file_y[file_y['Set']=="train"]['barcode'].tolist(), barcodes)
    elif "test" in file_y['Set'].tolist():  # if only test set is defined
        print("Setting test set from y-file")
        test_idx = indexInList(file_y[file_y['Set']=="test"]['barcode'].tolist(), barcodes)
    else:   # if nothing is defined
        print("Found Set column in Y-file but no mention of test set")

# Was the test set defined?
test_def = test_idx != []
# if the previous step didn't define everything, than randomly pick elements
if len(test_idx) + len(val_idx) + len(train_idx) != len(barcodes):
    # Number of test cells (this is relative to the full number of cells)
    nrTestCells = int(len(barcodes)*args.testSet)
    assert isinstance(args.maxCells, int)
    if nrTestCells > args.maxCells: nrTestCells = args.maxCells
    
    # Group cells by output
    test_groups = ["".join([str(i) for i in x]) for x in np.transpose(fullY.astype("int")).tolist()]
    test_groups_list={}
    for test_grp in set(test_groups):
        test_groups_list[test_grp] = []
    
    # Go through each group, identify indices for this group
    for idx, x in enumerate(test_groups):
        if not idx in train_idx: # if the test set is defined, do not include those barcodes in the further steps
            test_groups_list[x].append(idx)
    
    # Split indices into test and train set
    print(outputs)
    test_idx = [] if not test_def else test_idx # if the test_idx is already defined then we keep this
    val_idx = []
    train_idx = []
    for test_grp in sorted(set(test_groups)):
        # Draw a number from all indices for test and validation set, then draw from those for the validation set only
        # if the test set is predefined, then we still draw barcodes for it, but then do not use them
        test_and_val_group_N = int(((len(test_groups_list[test_grp]) + 0.0)/len(test_groups)) * nrTestCells * 2)
        test_and_val_barcodes=np.random.choice(test_groups_list[test_grp], test_and_val_group_N, replace=False).tolist()
        val_idx_x = np.random.choice(test_and_val_barcodes, int(test_and_val_group_N/2), replace=False).tolist()
        test_idx_x = sorted(set(test_and_val_barcodes) - set(val_idx_x)) if not test_def else [] # if the test indices are already defined, then we add nothing
        train_idx_x = sorted(set(test_groups_list[test_grp]) - set(val_idx_x) - set(test_idx_x))
        # print(and assertions)
        print(test_grp + ": test: " + str(len(test_idx_x)) + " val: " + str(len(val_idx_x)) + " train: " + str(len(train_idx_x)) + " of:  " + str(len(test_groups_list[test_grp])))
        assert len(test_idx_x) + len(val_idx_x) + len(train_idx_x) == len(test_groups_list[test_grp])
        assert len(set(test_idx_x) & set(val_idx_x)) == 0
        assert len(set(test_idx_x) & set(train_idx_x)) == 0
        assert len(set(train_idx_x) & set(val_idx_x)) == 0
        # Add up lists (across test_groups)
        test_idx = test_idx + test_idx_x
        val_idx = val_idx + val_idx_x
        train_idx = train_idx + train_idx_x

# Final assertions for the split
assert len(test_idx) + len(val_idx) + len(train_idx) == len(barcodes), "Error assigning test, training, and validation set"
assert len(set(test_idx) & set(val_idx)) == 0, "Error assigning test, training, and validation set"
assert len(set(test_idx) & set(train_idx)) == 0, "Error assigning test, training, and validation set"
assert len(set(train_idx) & set(val_idx)) == 0, "Error assigning test, training, and validation set"

# assign test and training set
y_train = fullY[:,train_idx]
y_val = fullY[:,val_idx]
y_test = fullY[:,test_idx]
open(os.path.join(outPath, "barcodes_test.txt"),"w").write("\n".join([barcodes[i] for i in test_idx]))
x_train = fullData[:,train_idx]
x_val = fullData[:,val_idx]
x_test = fullData[:,test_idx]

# print(result of draws)
print("Training Ys \t(total " + str(y_train.shape[1]) + "): \t" + "(== 1) - ".join(outputs) + "(== 1): \t" + " - ".join([str(y_train[i,:].sum()) for i in range(y_train.shape[0])]))
print("Validation Ys \t(total " + str(y_test.shape[1]) + "): \t" + "(== 1) - ".join(outputs) + "(== 1): \t" + " - ".join([str(y_test[i,:].sum()) for i in range(y_test.shape[0])]))
print("Testing Ys \t(total " + str(y_val.shape[1]) + "): \t" + "(== 1) - ".join(outputs) + "(== 1): \t" + " - ".join([str(y_val[i,:].sum()) for i in range(y_val.shape[0])]))


# Normalize to 0-1 for each row (input/gene) --------------------------------------------------------------------------------------------------------------------
if args.normalizedMatrix01:
    x_train = x_train.tocsr()
    x_val = x_val.tocsr()
    x_test = x_test.tocsr()
    # Go through all rows
    assert x_train.shape[0] == x_val.shape[0]
    assert x_train.shape[0] == x_test.shape[0] 
    for i in range(fullData.shape[0]):
        x_train_row=x_train.data[x_train.indptr[i]:x_train.indptr[i+1]]
        x_val_row=x_val.data[x_val.indptr[i]:x_val.indptr[i+1]]
        x_test_row=x_test.data[x_test.indptr[i]:x_test.indptr[i+1]]
        # if all cells in this row have a value (can be 0)
        if x_train_row.shape[0] == x_train.shape[1]:
            x_train_row_min=x_train_row.min()
            x_train_row=x_train_row-x_train_row_min
            x_val_row=x_val_row-x_train_row_min
            x_test_row=x_test_row-x_train_row_min
        
        if x_train_row.shape[0] > 0:    # if there are values in this row in the training data, devide by max of training data. 
            x_train_row_max = x_train_row.max()
            x_train_row=x_train_row/x_train_row_max
            x_val_row=x_val_row/x_train_row_max
            x_test_row=x_test_row/x_train_row_max
        else:
            if x_val_row.shape[0] > 0:   # Else but if there are values in validation, then normalize test data separately
                if x_val_row.shape[0] == x_val.shape[1]: # if all validation data have values, substract minimum (very unlikely but for completeness sake)
                    x_val_row=x_val_row-x_val_row.min()
                x_val_row=x_val_row/x_val_row.max()
            if x_test_row.shape[0] > 0:   # Else but if there are values in test, then normalize test data separately
                if x_test_row.shape[0] == x_test.shape[1]: # if all test data have values, substract minimum (very unlikely but for completeness sake)
                    x_test_row=x_test_row-x_test_row.min()
                x_test_row=x_test_row/x_test_row.max()
        
        x_train.data[x_train.indptr[i]:x_train.indptr[i+1]]=x_train_row
        x_val.data[x_val.indptr[i]:x_val.indptr[i+1]]=x_val_row
        x_test.data[x_test.indptr[i]:x_test.indptr[i+1]]=x_test_row


# Get weight matrix FUNCTION  --------------------------------------------------------------------------------------------------------------------
def weightMatrixFromYs(ys_input):
    matrix_groups = ["".join([str(i) for i in x]) for x in np.transpose(ys_input.astype("int")).tolist()]
    matrix_groups_weight={}
    # get factor to multiply by for each group
    for matrix_grp in sorted(set(matrix_groups)):
        matrix_groups_weight[matrix_grp] = (1.0/((matrix_groups.count(matrix_grp) + 0.0) / len(matrix_groups))) / len(set(matrix_groups))
    
    # create a numpy array with the same dimensions as the input
    return(np.array([[matrix_groups_weight[x] for x in matrix_groups],]* ys_input.shape[0]))

# Prepare TEST data --------------------------------------------------------------------------------------------------------------------
x_val = x_val.toarray()
y_val_weights = weightMatrixFromYs(y_val)
x_test = x_test.toarray()
y_test_weights = weightMatrixFromYs(y_test)

# SHUFFLE DATA --------------------------------------------------------------------------------------------------------------------
if args.shuffleGenes:
    rd.shuffle(genesList)
    print("Genes were shuffled")

if genesList[1] != genesListOrig[1]:
    print("Genes were indeed shuffled")

# FINISHED DATA SETUP --------------------------------------------------------------------------------------------------------------------
logMem("Setup data")



##############################
## MINIBATCH SIZE ############
##############################
# MINIBATCH SIZE ADJUSTMENT TO MAKE THEM EQUAL SIZE (otherwise the last minibatch is very small)
if args.minibatch == 0 or args.minibatch > x_train.shape[1]:
    args.minibatch = x_train.shape[1]
else:
    args.minibatch = int(x_train.shape[1]//(x_train.shape[1]//(args.minibatch)))



####################################
#### RANK NODES IN NETWORK RANKS ###
####################################
nodesRanks = []
edges_ranks = copy.deepcopy(edgelistFile)
edges_ranks_rows = edges_ranks.shape[0] + 1
i = 1
# Cycle through edgelist and add leaves to nodesRanks, thus building up nodesRanks from the bottom
while edges_ranks.shape[0] < edges_ranks_rows: # stop when no nodes (leaves) were added to the nodesRanks, edges_ranks_rows is edges_ranks.shape[0] at the last iteration
    edges_ranks_rows = edges_ranks.shape[0]
    leaves = sorted(set(edges_ranks['child'].tolist()) - set(edges_ranks['parent'].tolist())) # Nodes that are only children (not parents) are leaves
    if i > 1:
        nodesRanks.extend(leaves) # leaves are added to the nodes ranks
    else:
        assert set(leaves).issubset(set(genesListOrig)), "Some genes not found in the data!"
    edges_ranks = edges_ranks.loc[~edges_ranks['child'].isin(leaves)] # remove interactions of leaves
    i += 1

assert edges_ranks.shape[0] == 0, "Edgelist is circular!" # if there are edges remaining after the above
# add output nodes to the list
nodesRanks.extend(set(edgelistFile['parent'].tolist()) - set(edgelistFile['child'].tolist()))



###################################
#### MAP CONNECTIONS IN NETWORK ###
###################################
print("----------------\n...Mapping connections")
# edge list (for each node, get all the children (as string) in a list)
edges={}
for index, row in edgelistFile.iterrows():
    if(row['parent'] in edges.keys()):
        edges[row['parent']] = edges[row['parent']] + [row['child']]
    else:
        edges[row['parent']] = [row['child']]

# create nodeGeneMap and nodeNodeMap (for each node, has a list of all children genes / nodes in a list (as indices of genesList / nodesRanks))
# also map nodes to weights
nodeGeneMap = {}   # contains the gene from each node
nodeNodeMap = {}   # contains the nodes from each node
weightMap = {}     # the idx of the weights used (for each node a list of length genes + nodes)
iIt = 0
for iNode in nodesRanks:
    nodeGeneMap[iNode] = []
    nodeNodeMap[iNode] = []
    weightMap[iNode]=[]
    for eg in edges[iNode]:
        weightMap[iNode].append(iIt)
        iIt += 1
        if eg in genesList:
            nodeGeneMap[iNode].append(genesList.index(eg))
        elif eg in nodesRanks:
            nodeNodeMap[iNode].append(nodesRanks.index(eg))
        else:
            sys.exit(eg + " was not found in the nodes or the genes! Exiting")
            #print((eg + " was not found in the nodes or the genes! Exiting"))

# Weights
weightTotalLength = sum([len(weightMap[i]) for i in nodesRanks])

# Assertions for this set up
assert max([max(x) for x in weightMap.values()]) == weightTotalLength - 1
for iNode in nodesRanks:
    tarLength = len(nodeGeneMap[iNode]) + len(nodeNodeMap[iNode])
    assert tarLength == len(weightMap[iNode])
    assert tarLength == len(edges[iNode])
    assert len(nodeNodeMap[iNode]) == len(list(set(edges[iNode]) & set([nodesRanks[i] for i in nodeNodeMap[iNode]])))
    assert len(nodeGeneMap[iNode]) == len(list(set(edges[iNode]) & set([genesList[i] for i in nodeGeneMap[iNode]])))

# DONE
logMem("Setup Network")



#########################
#### TF NETWORK SETUP ###
#########################
print("----------------\n...Defining computational graph")

# Dropout Keep probability --------------------------------------------------------------------------------------------------------------------
dropoutKP_NODES = tf.cast(tf.placeholder_with_default(input=1.0, shape=[], name="dropoutKP_NODES"), tf.float64)
dropoutKP_GENES = tf.cast(tf.placeholder_with_default(input=1.0, shape=[], name="dropoutKP_GENES"), tf.float64)

# OUTPUT/Y --------------------------------------------------------------------------------------------------------------------
y_true = tf.placeholder(name="y_true", shape=[len(outputs),None],dtype=tf.float64)
y_weights = tf.placeholder(name="y_weights", shape=[len(outputs),None],dtype=tf.float64)

# INPUT/Genes --------------------------------------------------------------------------------------------------------------------
with tf.name_scope("Input_genes"):
    genesOrig = tf.placeholder(name="genes",shape=[x_train.shape[0], None],dtype=tf.float64)
    if doDropoutGenes:
        genes = tf.nn.dropout(x=genesOrig, keep_prob=dropoutKP_GENES)
    else:
        genes = genesOrig

# Nodes to approximate numerical gradients --------------------------------------------------------------------------------------------------------------------
numApproxPlus = tf.placeholder_with_default(input="", shape=[], name="NumericalApproximationUp_Node")
numApproxMinus = tf.placeholder_with_default(input="", shape=[], name="NumericalApproximationDown_Node")
numApproxVec = tf.reshape(tf.cast(tf.tile([args.numGradEpsilon], [tf.shape(genesOrig)[1]],), tf.float64), [1, tf.shape(genesOrig)[1]])

# Edge weights and intercept --------------------------------------------------------------------------------------------------------------------
with tf.name_scope("weights"):
    weights = tf.Variable(tf.random_normal([weightTotalLength,1], dtype=tf.float64, name="Random_weights"), name="Weights", dtype=tf.float64)
    tf.summary.scalar("mean", tf.reduce_mean(weights))
    tf.summary.histogram("histogram", weights)

with tf.name_scope("intercept"):
    interceptWeights = tf.Variable(tf.random_normal([len(nodesRanks)], dtype=tf.float64, name="Random_Intercepts"), name="Intercepts", dtype=tf.float64)
    tf.summary.scalar("mean", tf.reduce_mean(interceptWeights))
    tf.summary.histogram("histogram", interceptWeights)

# Hook up the nodes --------------------------------------------------------------------------------------------------------------------
nodes = {}
nodes_loss = {}
genes_unstacked = tf.unstack(genes)
for n in nodesRanks:
    #print(n)
    with tf.name_scope(n):
        # Z
        weightsX = tf.slice(weights, [weightMap[n][0],0],[len(weightMap[n]),1])
        features = tf.stack([genes_unstacked[x] for x in nodeGeneMap[n]] + [nodes[nodesRanks[nidx]] for nidx in nodeNodeMap[n]])
        z = tf.matmul(tf.transpose(weightsX), features) 
        z = z + tf.slice(interceptWeights, [nodesRanks.index(n)],[1])
        
        if not n in outputs:
            # Sigmoid, not on output nodes --> to be used in sigmoid_cross_entropy_with_logits function
            nodes[n] = tf.nn.sigmoid(z)
        else:
            nodes[n] = z
        
        # Numerical gradient approximation for this node?
        nodes[n] = tf.cond(tf.equal(numApproxPlus, tf.constant(n)), lambda: nodes[n] + numApproxVec, lambda: nodes[n])
        nodes[n] = tf.cond(tf.equal(numApproxMinus, tf.constant(n)), lambda: nodes[n] - numApproxVec, lambda: nodes[n])
        # unstack for later (parent nodes)
        nodes[n] = tf.unstack(nodes[n])[0]
        
        # Dropout, not on output nodes
        if doDropout and not n in outputs:
            print("dropout on " + n)
            # The next code tests if the parents of node n have more than one children. If this is not the case, then we do use dropout
            parents = edgelistFile[edgelistFile['child'].str.match("^" + n + "$")]['parent'].tolist()
            children = len(set(edgelistFile[edgelistFile['parent'].isin(parents)]["child"].tolist()))
            if children == 1:
                print("Skipping dropout for " + n)
            elif children == 2 and not args.disableDropoutAdjust:
                print(n + "'s parents have 2 children - dropout adjusted to 0.9 or " + str(args.dropOut))
                nodes[n] = tf.nn.dropout(x=nodes[n], keep_prob=tf.maximum(dropoutKP_NODES, 0.9))
            elif children == 3 and not args.disableDropoutAdjust:
                print(n + "'s parents have 3 children - dropout adjusted to 0.7 or " + str(args.dropOut))
                nodes[n] = tf.nn.dropout(x=nodes[n], keep_prob=tf.maximum(dropoutKP_NODES, 0.7))
            else:
                nodes[n] = tf.nn.dropout(x=nodes[n], keep_prob=dropoutKP_NODES)

# DONE with Network
logMem("Setup Tensorflow NW")

# Follow cost and other parameters
with tf.name_scope("regularizationCost"):
    weightsToPenalize = weights
    regularization = (args.lambd * tf.nn.l2_loss(weightsToPenalize)) / tf.cast(tf.shape(y_true)[1], tf.float64)
    tf.summary.scalar("value", regularization)

with tf.name_scope("Xentropy"):
    xentropy = tf.nn.sigmoid_cross_entropy_with_logits(logits=tf.stack([nodes[x] for x in outputs]), labels=y_true) * y_weights
    tf.summary.scalar("mean", tf.reduce_mean(xentropy))
    tf.summary.histogram("histogram", xentropy)

with tf.name_scope("Loss"):
    loss = tf.reduce_mean(xentropy) + regularization
    tf.summary.scalar("value", loss)

with tf.name_scope("y_hat"):
    y_hat = tf.nn.sigmoid(tf.stack([nodes[x] for x in outputs]))

with tf.name_scope("error"):
    error = tf.abs(y_true - y_hat) * y_weights
    tf.summary.scalar("mean", tf.reduce_mean(error))

with tf.name_scope("accuracy"):
    accuracy = tf.reduce_mean(tf.cast(tf.equal(tf.round(y_hat), y_true), tf.float64))
    tf.summary.scalar("mean", tf.reduce_mean(accuracy))

# Define training optimizer
with tf.name_scope("Optimizer"):
    if doADAM:
        optimizer = tf.train.AdamOptimizer(learning_rate=args.alpha)
    else:
        optimizer = tf.train.MomentumOptimizer(learning_rate=args.alpha, momentum=args.momentum)
    
    train = optimizer.minimize(loss)

# DONE
logMem("Setup Tensorflow")



###################
#### INITIALIZE ###
###################
print("----------------\n...Initializing")
init = tf.global_variables_initializer()
sess = tf.Session(config=tf.ConfigProto(inter_op_parallelism_threads=args.threads,intra_op_parallelism_threads=args.threads))
sess.run(init)
logMem("Initiated Tensorflow")



#############################
#### PREPARE OUTPUT FILES ###
#############################
print("----------------\n...Preparing outputs")
sep=","

# Saver (saves models during training)
saver = tf.train.Saver()

# Summary writers (Tensorboard)
merged = tf.summary.merge_all()
if args.tfWrite:
    train_writer = tf.summary.FileWriter(outPath + '/train', sess.graph)
    test_writer = tf.summary.FileWriter(outPath + '/test', sess.graph)

# Cost output (Track progress)
costFile = open(os.path.join(outPath, "tf_cost.csv"),"w")
costFile.write("iteration,train,validation,error,accuracy,breakCounter" + "\n")
costFile.close()

# Settings
settingsFile = open(os.path.join(outPath, "tf_settings.csv"),"w")
settingsFile.write("alpha" + sep + str(args.alpha) + "\n")
settingsFile.write("lambda" + sep + str(args.lambd) + "\n")
settingsFile.write("DropoutKP" + sep + str(args.dropOut) + "\n")
settingsFile.write("DropoutKPGenes" + sep + str(args.dropOutGenes) + "\n")
settingsFile.write("momentum" + sep + str(args.momentum) + "\n")
settingsFile.write("minibatch" + sep + str(args.minibatch) + "\n")
settingsFile.write("python" + sep + sys.version + "\n")
settingsFile.write("Seed" + sep + str(args.randomSeed) + "\n")
settingsFile.close()



#################
#### TRAINING ###
#################
print("----------------\n...Starting training")
# training iterations
trainCost = []
testErr = []
testErrCurrentMin = 1
breakCounter = 0
start = time.time()
logMem("Prepared Training")
for i in [xx + 1 for xx in range(args.iterations)]:
    
    # DEFINE MINIBATCH    
    idxs = list(range(x_train.shape[1]))
    np.random.shuffle(idxs)
    minibatch_list = [idxs[j:min(j+args.minibatch, len(idxs)-1)] for j in list(range(0, x_train.shape[1], args.minibatch))]
    # generates list of lists, each containing the indices used
    
    # TRAIN BY MINIBATCH
    for i_batch in range(len(minibatch_list)):
        # First train for this epoch
        batch_idx=minibatch_list[i_batch]
        x_batch=x_train[:,batch_idx].toarray()
        y_batch=y_train[:,batch_idx]
        y_batch_weights = weightMatrixFromYs(y_batch)
        if i > 1: # The first iteration reflects initiation parameters (no training done)
            sess.run(train, {genesOrig:x_batch, y_true:y_batch, dropoutKP_NODES: args.dropOut, dropoutKP_GENES: args.dropOutGenes, y_weights:y_batch_weights})
    
    if i < 10 or i % 10 == 0:
        print("\nTraining epoch: " + str(i))
        
        # Loss, error, etc on training and validation set
        lossRes, trainClassProb, trainWriteContent = sess.run([loss, y_hat, merged], {genesOrig:x_batch, y_true:y_batch, dropoutKP_NODES: args.dropOut, dropoutKP_GENES: args.dropOutGenes, y_weights:y_batch_weights})
        lossTestRes, testWriteContent, testerror, testAccuracy = sess.run([loss, merged, error,accuracy], {genesOrig:x_val, y_true:y_val, y_weights:y_val_weights})
        testErrRun = testerror.mean()
        
        # print(and write progress)
        print("Mean loss: " + "%.4f" % lossRes + " Validation loss: " + "%.4f" % lossTestRes + " Validation error: " + str(testErrRun))
        costFile = open(os.path.join(outPath, "tf_cost.csv"),"a")
        costFile.write(str(i) + sep + str(lossRes) + sep + str(lossTestRes) + sep + str(testErrRun) + sep + str(testAccuracy.mean()) + sep + str(breakCounter) + "\n")
        costFile.close()
        
        # Tensorboard output
        if args.tfWrite:
            train_writer.add_summary(trainWriteContent, i)
            test_writer.add_summary(testWriteContent, i)
        
        # Early stopping
        if i > 10 and not args.disableInterrupt:
            if trainCost[-1] - lossRes < 0.00001: # if the general error goes up or stays equal
                breakCounter += 1
            if any(testErrPrevious + 0.05 < testErrRun for testErrPrevious in testErr): # if test error goes up
                breakCounter += 1
        
        print("Break Counter = " + str(breakCounter))
        
        if breakCounter > args.maxBreakCount:
            break
        
        # Save model
        if i > 2 and (testErrRun < 0.2 and testErrCurrentMin * args.minImprovement > testErrRun):
            breakCounter = 0
            saver.save(sess, os.path.join(outPath,'best-model'))
            testErrCurrentMin = testErrRun
            print("Model saved at iteration " + str(i) + " at Test error: " + str(testErrRun))
        
        trainCost.append(lossRes)
        testErr.append(testErrRun)
    

# DONE
print("----------------\n...Training done:")
print(start - time.time())
logMem("Trained Tensorflow")



#########################################
#### OUTPUT WEIGHTS, yHat, yTrue, ... ###
#########################################
# Load last model
if os.path.exists(os.path.join(outPath,'best-model.meta')):
    saver.restore(sess, os.path.join(outPath,'best-model'))
    print("Model restored from " + outPath)

logMem("Loaded Best Model")

# Write out Test Error with loaded model
valErrNumGrad = sess.run(error, {genesOrig:x_test, y_true:y_test, y_weights:y_test_weights}).mean()
print("Test Error before Numerical Gradient = " + str(valErrNumGrad))
open(os.path.join(outPath, "tf_NumGradTestError.txt"),"w").write("NumGradTestError" + "\n" + str(valErrNumGrad))


# Write y hat to file
val_y_hat = sess.run(y_hat, {genesOrig:x_test, y_true:y_test})
print("----------------\n...Writing y hat to file")
val_y_hat_file = open(os.path.join(outPath, "tf_yHat_test.csv"),"w")
val_y_hat_file.write(",".join(outputs) + "\n")
val_y_hat_file.write(pd.DataFrame(np.transpose(val_y_hat)).to_csv(header=False, index=False))
val_y_hat_file.close()

# Write y true to file
val_y_hat_file = open(os.path.join(outPath, "tf_yTrue_test.csv"),"w")
val_y_hat_file.write(",".join(outputs) + "\n")
val_y_hat_file.write(pd.DataFrame(np.transpose(y_test.astype("int"))).to_csv(header=False, index=False))
val_y_hat_file.close()

logMem("Tensorflow results accuracy")

# Weight output
sep = ","
weightRes, interceptRes = sess.run([weights, interceptWeights])
print("weights: " + "%.4f" % np.abs(weightRes).mean() + " intercepts: " + "%.4f" % np.abs(interceptRes).mean())
weightFile = open(os.path.join(outPath, "tf_weights.csv"),"w")
weightFile.write("parent,child,weight" + "\n")
for n in nodesRanks:
    i = 0
    for targetGene in nodeGeneMap[n]:
        weightFile.write(n + sep + genesList[targetGene] + sep + str(weightRes[weightMap[n][i]].tolist()[0]) + "\n")
        i += 1
    for targetProtein in nodeNodeMap[n]:
        weightFile.write(n + sep + nodesRanks[targetProtein] + sep + str(weightRes[weightMap[n][i]].tolist()[0]) + "\n")
        i += 1
    assert i == len(weightMap[n])
    weightFile.write(n + sep + "Intercept" + sep + str(interceptRes[nodesRanks.index(n)]) + "\n")

weightFile.close()



###########################################
#### NUMERICAL GRADIENT NODE IMPORTANCE ###
###########################################
if not args.disableNumGrad:
    print("----------------\n...Starting Numerical Gradient estimation:")
    print(start - time.time())

    #numGrad_yWeights = np.array([[1.0 for x in range(y_test.shape[1])],]* y_test.shape[0]) # REMOVE IF HAS NOT FAILED
    # For each node, get y_hat with pos and neg perturbation (numApproxPlus/numApproxMinus)
    # calculate mean across 
    numAgg = 0
    nodesNumGrad = [x for x in nodesRanks if not x in outputs]
    for n in nodesNumGrad:
        print("Calculating numerical gradient for: " + n)
        logMem("Tensorflow numgrad " + n)
        res = (sess.run(y_hat, {genesOrig:x_test, y_true:y_test, numApproxPlus: n}) - sess.run(y_hat, {genesOrig:x_test, y_true:y_test, numApproxMinus: n}))/(2*args.numGradEpsilon)
        res = res.mean(1)
        if numAgg.__class__ != np.ndarray:
            numAgg = res
        else:
            numAgg = np.column_stack((numAgg, res))

    # Write means to file
    sep = ","
    numGradFile = open(os.path.join(outPath, "tf_NumGradMeans.csv"),"w")
    numGradFile.write("Node" + sep + sep.join(outputs) + "\n")
    for i,n in enumerate(nodesNumGrad):
        numGradFile.write(n + sep + sep.join([str(x) for x in numAgg[:,i].tolist()]) + "\n")

    numGradFile.close()

    # done with NUMGRAD
    logMem("Tensorflow numgrad done")



##################
#### FINISH UP ###
##################
print("----------------\n...Time passed:")
print(start - time.time())

sess.close()

print("\n\n\t\t\tKPNN TRAINING COMPLETED SUCCESSFULLY\n\n")