# xor-haplotyping

This package contains the MATLAB implementation of the proposed xor-haplotyping (XHSD) algorithm. 

Please see the script textXHSD.m to test XHSD algorithm on CFTR data with the desired parameters.

The main program is sparsehaplotypeSPL.m, which demands as input the genotype data. The genotypes should be located in the columns of the input matrix.

INSTRUCTIONS for saprsehaplotypeSPL.m:

Input Format: L x N matrix corresponding to N individuals with L SNPs that are typed as follows.<br/> 
0: homozygous-common 00<br/>
2: homozygous-mutant 11<br/>
1: heterozygous 01 or 10<br/>
4: homozygous-hidden (XOR-site) 00 or 11<br/>
6: missing data

Parameters:<br/>
W = maximum block size in PL method<br/>
interest = which individuals to infer haplotypes. Use [1, 2, .., N] as default.<br/> 
Lmiss= L x N binary matrix where 1 means a missing site.<br/>
flagBlock = 1:PL method (default), 2:fixed-length partitioning, 0:no-partition<br/>

Output Format: L x N x 2 matrix.<br/>
Haplotypes are located in the corresponding columns in dimension 1 (L,N,1) and dimension 2 (L,N,2).<br/>
￼￼￼￼￼

Please cite:<br/>
Abdulkadir Elmas, Guido H Jajamovich and Xiaodong Wang<br/>
Maximum parsimony xor haplotyping by sparse dictionary selection<br/>
BMC Genomics 2013 14:645<br/>
https://doi.org/10.1186/1471-2164-14-645
