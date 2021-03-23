# Two-Step Algorithm

Matlab functions to run the Two-Step algorithm.

A dataset X of non-Gaussian variables is required as input, together with a positive value for the penalization parameter, lambda.

Two-Step outputs the causal coefficients matrix B, from X = BX + E. 

In B, the causal direction goes from column to row, such that a matrix entry Bij, implies Xj --> Xi

two_step_CD.m  runs Two-Step with Adaptative Lasso as a first step. 

two_step_CD_mask.m  should be used if the adjacency matrix of the first step was computed with some other algorithm.
