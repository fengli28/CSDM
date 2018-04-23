The CSDM method is proposed to detect cancer specific driver modules of a certain cancer based on other cancer types
by using the mutation data.
The MATLAB code of the CSDM method can be downloaded from 

The datasets used in CSDM are provided in the folder "\CSDM\0 Mutation matrix" including 11 cancer types.
The mutation matrix is a binary matrix and has n rows (samples) and m columns (samples).
For example:
	TCGA-FD-A3NA	TCGA-FD-A3N5	TCGA-FD-A3N6	TCGA-DK-A1AF	TCGA-BT-A0S7	TCGA-DK-A3IV	
A1BG	0	1	0	0	0	0	
A2M	0	0	0	0	1	0	


The operation method is as follows:
1. Calculate the specific coverage and mutual exclusivity for each cancer type.
Under the folder "\CSDM\1 Calculate CO ME", open the file "MainCOME.m", and give the following files.
Two inputs:
---AllCancer.txt (the mutation matrix for all cance types)
---BLCA.txt (the mutation matrix for cancer k)

Run:
MainCOME

Output:
---CO_BLCA.mat (the specific coverage matrix of cancer k)
---ME_BLCA.mat (the mutual exclusivity matrix of cancer k)

After the program "MainCOME" runs out, we get two matrices, coverage matrix CO_* and mutual exclusivity matrix ME_*
for a certain cancer.
This step should run on each cancer type to obtain coverage matrix and mutual exclusivity matrix for each cancer.

2. Sort gene pairs across all cancer types. 
1) Specific coverage
Under the folder "\CSDM\2 Sort gene pairs", open the file "CompareAllCO.m", and give the following files.
K inputs:
---CO_BLCA.mat (the specific coverage matrix of cancer 1)
---CO_BRCA.mat (the specific coverage matrix of cancer 2)
---...
---CO_UCEC.mat (the specific coverage matrix of cancer K)

Run
CompareAllCO

K Outputs:
---TopCO_BLCA (a binary matrix represents the gene pairs with the largest specific coverage in cancer 1 among all cancers)
---TopCO_BRCA (a binary matrix represents the gene pairs with the largest specific coverage in cancer 2 among all cancers)
...
---TopCO_UCEC (a binary matrix represents the gene pairs with the largest specific coverage in cancer K among all cancers)

When the program "CompareAllCO" runs out, we get 11 coverage rank matrices for all cancer types, denoted by TopCO_*. 
The coverage rank matrix for cancer k is a binary matrix, and the element is one if the corresponding gene pair has the 
largest coverage across all cancer types. 

2) Mutual exclusivity
Under the folder "\CSDM\2 Sort gene pairs", open the file "CompareAllME.m", and give the following files.
K inputs:
---ME_BLCA.mat (the mutual exclusivity matrix of cancer 1)
---ME_BRCA.mat (the mutual exclusivity matrix of cancer 2)
---...
---ME_UCEC.mat (the mutual exclusivity matrix of cancer K)

Run
CompareAllME

K Outputs:
---TopME_BLCA (a binary matrix represents the gene pairs with the largest mutual exclusivity in cancer 1 among all cancers)
---TopME_BRCA (a binary matrix represents the gene pairs with the largest mutual exclusivity in cancer 2 among all cancers)
...
---TopME_UCEC (a binary matrix represents the gene pairs with the largest mutual exclusivity in cancer K among all cancers)

When the program "CompareAllME" runs out, we get 11 mutual exclusivity rank matrices for all cancer types, denoted by TopME_*. 
The mutual exclusivity rank matrix for cancer k is a binary matrix, and the element is one if the corresponding 
gene pair has the largest mutual exclusivity across all cancer types. 

3. Specific network construction and specific modules detection.
Under the folder "\CSDM\3 Specific driver module detection", open the file "main.m", and give the following files.
Six inputs:
---AllCancer.txt (the mutation matrix for all cance types)
---BLCA.txt (the mutation matrix for cancer k)
---TopCO_BLCA.mat (a binary matrix denotes the gene pairs with the largest coverage in cancer k)
---TopME_BLCA.mat (a binary matrix denotes the gene pairs with the largest mutual exclusivity in cancer k)
---CO_BLCA.mat (the specific coverage matrix of cancer k)
---ME_BLCA.mat (the mutual exclusivity matrix of cancer k)

ran

main

Output:
Modules of cancer k.txt

When the program "main" runs out, we get a file "Modules of cancer k.txt", which records the specific modules of cancer k. 
For example:
Number of genes in module	Specific coverage	Internal coverage	External coverage	Genes in module
8	0.521	0.609	0.445	AMIGO3	CDKAL1	CDKN1A	ELF3	GGA3	NPIP	SCNN1A	SPECC1L	
.
.
.



