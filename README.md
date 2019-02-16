# Smith-Waterman 

## Overview
This is a basic implementation of Smith-Waterman local alignment algorithm for MCDB 452: Biomedical Data Science.

## Install
Quick method to install:
    
    sudo pip install git+git://github.com/nszheng/Smith-Waterman.git
    
## Usage
    python sw.py -i <input file> -s <score file> -o <open gap weight> -e <extend gap weight>
    
## Example Code
    python sw.py -i sample-input.txt -s blosum62.txt -o -2 -e -1

## Sample Input
Input .txt file has sequence 1 on line 1 and sequence 2 on line 2:
```
FDKFKHLK
KLFPKFAGIAHGDL
```

## Sample Output
Sample output produced using the blosum62.txt score matrix. It displays the input sequences, the score matrix, and the best local alignment with the alignment score. 

Additional example input and output with longer sequences are contained in the Examples folder. 

```
----------------------------------------------------
|	Sequences                                   |
----------------------------------------------------

sequence 1:
FDKFKHLK
sequence 2:
KLFPKFAGIAHGDL

----------------------------------------------------
|	Score Matrix                                |
----------------------------------------------------

	F	D	K	F	K	H	L	K	
	0	0	0	5	3	5	3	2	5
K	0	0	0	5	3	5	3	2	5
L	0	0	0	3	5	3	2	7	5
F	0	6	4	3	9	7	6	5	4
P	0	4	5	3	7	8	6	5	4
K	0	3	3	10	8	12	10	9	10
F	0	6	4	8	16	14	13	12	11
A	0	4	4	7	14	15	13	12	11
G	0	3	3	6	13	13	13	11	10
I	0	2	1	5	12	12	11	15	13
A	0	1	0	4	11	11	10	13	14
H	0	0	0	3	10	10	19	17	16
G	0	0	0	2	9	9	17	15	15
D	0	0	6	4	8	8	16	14	14
L	0	0	4	4	7	7	15	20	18


----------------------------------------------------
|	Best Local Alignment                         |
----------------------------------------------------

Alignment Score: 20

  (FDKF---KH--L)K
   | ||    |  |
KL(FPKFAGIAHGDL)
```
