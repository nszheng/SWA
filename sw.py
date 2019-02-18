#!/usr/bin/python
__author__ = "Neil Zheng"
__email__ = "neil.zheng@yale.edu"
__copyright__ = "Copyright 2019"
__license__ = "GPL"
__version__ = "1.0.0"
### Usage: python hw1.py -i <input file> -s <score file>
### Example: python hw1.py -i input.txt -s blosum62.txt
### Note: Smith-Waterman Algorithm
### Scripting must be done from scratch, without the use of any pre-existing packages.
### Python standard library (I/O) and numpy are allowed.
import argparse
import numpy as np
### This is one way to read in arguments in Python.


### We need to read input file and score file.
parser = argparse.ArgumentParser(description='Smith-Waterman Algorithm')
parser.add_argument('-i', '--input', help='input file', required=True)
parser.add_argument('-s', '--score', help='score file', required=True)
parser.add_argument('-o', '--opengap', help='open gap', required=False, default=-2)
parser.add_argument('-e', '--extgap', help='extension gap', required=False, default=-1)
args = parser.parse_args()


### Implement your Smith-Waterman Algorithm
def runSW(inputFile, scoreFile, openGap, extGap):
    # Opening input file and storing sequence as strings seq1 and seq2
    f = open(inputFile, "r")
    seq1 = f.readline()
    seq2 = f.readline()
    f.close()
    
    # Opening score file, storing column names as as a list simList 
    # and the similarity matrix as numpy matrix simMat
    f = open(scoreFile, "r")
    simList = f.readline().split()
    ncols = len(simList)
    f.close()
    # We remove column names and row names for the similarity matrix
    simMat = np.loadtxt(scoreFile, skiprows = 1, usecols = range(1,ncols))
    
    
    # Initializing a numpy array H to hold alignment matrix scores with 
    # additional dummy row and column of initial zeros for algorithm
    H = np.zeros((len(seq2), len(seq1)))
    
    # Initializing value and position of best alignment score
    (ibest,jbest) = (0, 0)
    bestscore = 0
    
    # Initializing arrays to hold the prior position of sequence
    priori = np.zeros((len(seq2), len(seq1)))
    priorj = np.zeros((len(seq2), len(seq1)))
    
    ###### Nested for loop to fill out alignMat ######
    for i in range(1, len(seq2)):
        for j in range(1, len(seq1)):
            # Calculating the value for H[i,j] if seq1[j-1] and seq2[i-1] similar
            sim = H[i-1,j-1] + simMat[simList.index(seq1[j-1]), simList.index(seq2[i-1])]
            
            # Calculating value of H[i,j] if deletion of length k ends at seq1[j-1]
            # We compute the max over the prior values of H in the row plus 
            # a opening gap penalty and multipled by a extension gap penalty.
            # We also store its position for later
            delseq1_pos = np.argmax(H[range(i),j] + 
                               extGap*np.flip(range(i), axis = 0) + 
                               openGap)
            delseq1 = H[delseq1_pos, j] + extGap*(i-1-delseq1_pos) + openGap
            
            # Calculating value of H[i,j] if deletion of length k ends at seq2[i-1]
            delseq2_pos = np.argmax(H[i,range(j)] + 
                               extGap*np.flip(range(j), axis = 0) + 
                               openGap)
            delseq2 = H[i,delseq2_pos] + extGap*(j-1-delseq2_pos) + openGap
            
            # Storing final value of H[i,j]
            H[i,j] = np.max([sim,delseq1,delseq2,0])
            
            
            # Updating prior sequence location
            # Finding index of max value
            maxInd = np.argmax([sim,delseq1,delseq2,0])
            # If sim was max, prior index is just i-1, j-1
            if maxInd == 0:
                priori[i,j] = i-1
                priorj[i,j] = j-1
            # If delseq1 was max, prior index is i, j-k
            elif maxInd == 1:
                priori[i,j] = delseq1_pos
                priorj[i,j] = j 
            # If delseq2 was max, prior index is i-k, j
            elif maxInd == 2:
                 priori[i,j] = i
                 priorj[i,j] = delseq2_pos
            # If 0 was max, no prior index
            else:
                priori[i,j] = i
                priorj[i,j] = j
            
            
        
            # Updating best score and its position
            if H[i,j] > bestscore:
                bestscore = H[i,j]
                (ibest,jbest) = (i,j)
                    
    
    ###### Determining best sequence via traceback ######
    # Initializing the first (i,j) to (ibest,jbest) for traceback 
    # and storing first aligned elements of sequence
    (i,j) = (ibest,jbest)
    alignseq1 = [seq1[j-1]]
    alignseq2 = [seq2[i-1]]         
    (inext,jnext) = (int(priori[i,j]), int(priorj[i,j]))
    
    ### While loop to implement traceback ###
    while H[inext,jnext] != 0:
        # Finding the (i,j) position and next (i,j) position
        (i,j) = (int(priori[i,j]), int(priorj[i,j]))
        (inext,jnext) = (int(priori[i,j]), int(priorj[i,j]))
        
        # If the jnext = j, then we have a deletion in sequence 1 and
        # we extend alignseq1 with '-' (i-inext) number of times and
        # continue appending (i-inext) times to alignseq2
        # We do similarly if inext == i.
        # Else, we append the next sequence element as specified by (i,j)
        if jnext == j:
            alignseq1.extend(["-"] * (i-inext))
            for k in range(0,i - inext):
                alignseq2.append(seq2[i-k-1])
        elif inext == i:
            for k in range(0,j - jnext):
                alignseq1.append(seq1[j-k-1])
            alignseq2.extend(["-"] * (j-jnext))
        else:
            alignseq1.append(seq1[j-1])
            alignseq2.append(seq2[i-1])
            
    # We will use inext and jnext to define the bounds of the non-aligned sequence
    
    # Reversing aligned sequence to be proper order
    alignseq1.reverse()
    alignseq2.reverse()
            
    
    ###### Printing output file ######
    f = open("output.txt", "w")
    
    ### Printing sequences ###
    f.write("----------------------------------------------------\n" + 
            "|	Sequences                                   |\n" +
            "----------------------------------------------------\n\n")
    f.write("sequence 1:\n" + seq1)
    f.write("sequence 2:\n" + seq2 + "\n")
    
    ### Printing score matrix ###
    f.write("----------------------------------------------------\n" +
            "|	Score Matrix                                |\n" +
            "----------------------------------------------------\n\n")
    mat = np.matrix(H) # Converting to matrix for easier printing
    f.write('\t' + '\t'.join(list(seq1))) # Writing sequence 1 for colnames
    
    (nrow, ncol) = np.shape(mat) # Getting dimensions of score matrix for printing
    
    # Printing first row of dummmy zeros
    f.write("\t") 
    np.savetxt(f, mat[1,:], fmt='%d', delimiter = "\t") 
    
    for i in range(1,nrow): # For Loop to print all rows in score matrix
        f.write(seq2[i-1] + '\t')
        np.savetxt(f, mat[i,:], fmt='%d', delimiter = "\t")
        
    
    ### Printing best local alignment ###
    f.write("\n\n----------------------------------------------------\n" + 
            "|	Best Local Alignment                         |\n" + 
            "----------------------------------------------------\n\n")
    f.write("Alignment Score: \t" + str(int(bestscore)) + "\n")
    
    
    # Determining which elements align perfectly for printing
    aligntrue = []
    for i in range(len(alignseq1)):
        if(alignseq1[i] == alignseq2[i]):
            aligntrue.append("|")
        else:
            aligntrue.append(" ")
    
    # Using inext and jnext from earlier, we can determine the bounds of the
    # non-aligned parts of the sequence and print them first
    f.write(" "*max(inext-jnext,0) + seq1[0:jnext] + 
            "(" + "".join(alignseq1) + ")" + seq1[jbest:len(seq1)])
    f.write (" "*max(inext+1, jnext+1) + "".join(aligntrue) + "\n")
    f.write(" "*max(jnext-inext,0) + seq2[0:inext] + 
            "(" + "".join(alignseq2) + ")" + seq2[ibest:len(seq2)]) 
    f.close()


    ### Print input and score file names. You can comment these out.
    print ("input file : %s" % inputFile)
    print ("score file : %s" % scoreFile)
    print ("open gap penalty : %s" % openGap)
    print ("extension gap penalty : %s" % extGap)
    
############################################## end runSW   


### Run your Smith-Waterman Algorithm
runSW(args.input, args.score, args.opengap, args.extgap)

