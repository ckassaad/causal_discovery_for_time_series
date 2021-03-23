README for Simulated financial time series data as described in Kleinberg, S (2013). Causality, Probability and Time. Cambridge University Press.

Version 10/16/2012

This data is made available under a Creative Commons Attribution-NonCommercial 3.0 Unported license (http://creativecommons.org/licenses/by-nc/3.0/).
-----------------------------------
INCLUDED FILES:

Directory returns/
manyinputs30007000_returns.csv
manyinputs800012000_returns.csv
nocause_returns30007000.csv
nocause_returns800012000.csv
random-rels_20_1A_returns30007000.csv
random-rels_20_1A_returns800012000.csv
random-rels_20_1B_returns30007000.csv
random-rels_20_1B_returns800012000.csv
random-rels_20_1C_returns30007000.csv
random-rels_20_1C_returns800012000.csv
random-rels_20_1D_returns30007000.csv
random-rels_20_1D_returns800012000.csv
random-rels_20_1E_returns30007000.csv
random-rels_20_1E_returns800012000.csv
random-rels_20_1_3A_returns30007000.csv
random-rels_20_1_3A_returns800012000.csv
random-rels_20_1_3B_returns30007000.csv
random-rels_20_1_3B_returns800012000.csv
random-rels_40_1A_returns30007000.csv
random-rels_40_1A_returns800012000.csv
random-rels_40_1B_returns30007000.csv
random-rels_40_1B_returns800012000.csv
random-rels_40_1_3A_returns30007000.csv
random-rels_40_1_3A_returns800012000.csv
random-rels_40_1_3B_returns30007000.csv
random-rels_40_1_3B_returns800012000.csv

Directory relationships/
manyinputs.csv			
random-rels_20_1_3.csv
random-rels_20_1A.csv
random-rels_20_1B.csv	
random-rels_20_1C.csv	
random-rels_20_1D.csv
random-rels_20_1E.csv		
random-rels_40_1.csv
random-rels_40_1_3.csv

figures/
manyinput.eps
manyinput.pdf
random20_1-3.eps
random20_1-3.pdf
random20_1A.eps
random20_1A.pdf
random20_1B.eps
random20_1B.pdf
random20_1C.eps
random20_1C.pdf
random20_1D.eps
random20_1D.pdf
random20_1E.eps
random20_1E.pdf
random40_1.eps
random40_1.pdf
random40_1_3.eps
random40_1_3.pdf

Files in /returns are the generated return time series. Files in /relationships are the actual embedded causal relationships for each return time series. Files in /figures are graphical representations of the sets of relationships.

FILE NAME CONVENTIONS:
random-rels refer to datasets where the relationships between portfolios was generated randomly.

random_rels_X_Y[_Z] means the dataset contains X relationships where each has a lag of at least Y time units between cause and effect. If Z is included the lags range from Y to Z (inclusive), otherwise all lags are exactly Y. 

30007000 and 800012000 refer to the two different time periods (created by taking days 3K-7K and 8K-12K from the original data).

The file with the true relationships embedded in a dataset is denoted by the same prefix as the datafile before the underscore. Note that the relationships are the same for both time periods with a given structure, thus there are 9 structure files (since one structure has no causes).

RELATIONSHIP FILES:
Format is rows of cause, effect, lag (no header), where this means that in each row the portfolio from the first column influences the portfolio in the second column at the time lag (in days) of the lag in the third column.

DATA FILES:
Each data file is a csv file with 25 columns and 4000 rows. The ith column denotes observations of the ith stock, with the jth row being the jth day in the observation sequence.

FIGURES:
Figure files depicting each set of embedded relationships are provided in both eps and pdf formats. Only portfolios involved in causal relationships are included in the graph (thus all 25 portfolios may not be depicted in all graphs). Solid edges indicate a lag of one day, double edges mean two days and dashed edges mean a lag of 3 days.


-----------------------------------
Data generation procedure

The data were generated using a factor model, with causality embedded through the idiosyncratic terms. 

The return of a portfolio i at time t is given by:
ri,t= sum_j beta_ij f_it' + epsilon_i,t

This is a sum over the three factors along with their weightings for each stock. Here t'=t. The epsilons are the randomly generated idiosyncratic terms. When a portfolio k is a cause of portfolio i at time lag t, then this is the sum of i's epsilon along with all such epsilon_k,t-l.

Factors were constructed using the Fama-French daily factors 1963-2011 along with the 5x5 size/book-to-market portfolios. 
