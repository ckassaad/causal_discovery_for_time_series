# causal_discovery_for_time_series
Package to test causal discovery algorithm on simulated and real data 


THIS SOURCE CODE IS SUPPLIED AS IS WITHOUT WAR
RANTY OF ANY KIND	 AND ITS AUTHOR AND THE JOURNAL OF
ARTIFICIAL INTELLIGENCE RESEARCH JAIR AND JAIRS PUBLISHERS AND DISTRIBUTORS	 DISCLAIM ANY AND ALL WARRANTIES	
INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE	 
AND ANY WARRANTIES OR NON INFRINGEMENT THE USER
ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS
SOURCE CODE	 AND NEITHER THE AUTHOR NOR JAIR	 NOR JAIRS
PUBLISHERS AND DISTRIBUTORS	 WILL BE LIABLE FOR DAM
AGES OF ANY KIND RESULTING FROM ITS USE Without limiting
the generality of the foregoing	 neither the author	 nor JAIR	 nor JAIR's 
publishers and distributors	 warrant that the Source Code will be errorfree	
will operate without interruption	 or will meet the needs of the user


## Methods

* NBCB: https://link.springer.com/chapter/10.1007/978-3-030-86486-6_28
* PCTMI: 
* TiMINO: https://proceedings.neurips.cc/paper/2013/file/47d1e990583c9c67424d369f3414728e-Paper.pdf
* VarLiNGAM: https://www.jmlr.org/papers/volume11/hyvarinen10a/hyvarinen10a.pdf
* PCMCICMIknn: http://proceedings.mlr.press/v124/runge20a/runge20a.pdf
* PCMCIParCorr: http://proceedings.mlr.press/v124/runge20a/runge20a.pdf
* oCSE: https://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.474.6986&rep=rep1&type=pdf
* tsFCI: 
* DyNoTears: https://arxiv.org/pdf/2002.00498.pdf
* GrangerMV: https://pubmed.ncbi.nlm.nih.gov/20481753/
* GrangerPW
* TCDF: 

Some algorithms are imported from other langauges such as R and Java


## Test

### To test algorithms on simulated data run:
python3 test_simulated_data.py method structure n_samples num_processor verbose

* method: causal dicovery algorithms, choose from [GrangerPW, GrangerMV, TCDF, PCMCICMIknn, PCMCIParCorr, oCSE, PCTMI, tsFCI, VarLiNGAM, TiMINO, Dynotears]
* structure: causal structure, choose from [fork, v_structure, diamond, 7ts2h]
* n_samples: number of timestamps
* num_processor: number of processors 

Example: python3 test_fmri.py "NBCB" "fork" 1000 1 1

### To test algorithms on fmri data run:
python3 test_simulated_data.py method num_processor verbose

Example: python3 test_fmri.py "NBCB" 1 1
