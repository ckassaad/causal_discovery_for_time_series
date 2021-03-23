# KernelGrangerCausality
Code to evaluate nonlinear Granger causality using the kernel trick to reduce complexity
This set of functions implements the Kernel Granger Causality, allowing to detect nonlinear dynamical interactions, as described in

Kernel Granger causality and the analysis of dynamical networks
D. Marinazzo, M.Pellicoro and S. Stramaglia
Physical Review E, 77, 056215 (2008)
http://journals.aps.org/pre/abstract/10.1103/PhysRevE.77.056215

https://arxiv.org/abs/0803.2881

Kernel method for nonlinear Granger causality
D. Marinazzo, M.Pellicoro and S. Stramaglia
Phys. Rev. Lett. Vol. 100 pag. 144103 (2008)
http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.144103

https://arxiv.org/pdf/0711.2643.pdf

The folder contains two toolboxes

Kernel causality last

LOO_crossvalidation for choosing the model order

unzip these two toolboxes with the other MATLAB toolboxes and add them to your path. 

You will find among others the scripts: 

test_KGC_KCL: run kernel causality on a simulated network of 5 nodes

test_biv_2n: test KGC on two coupled maps

test_var_eps_1: modulation in time of the coupling parameter for the two coupled maps (contains calls to Granger causality and correlation with sliding window)

Please do not hesitate to contact us for suggestions and remarks http://ost.io/@danielemarinazzo

DISCLAIMER OF WARRANTIES AND LIMITATION OF LIABILITY The code is supplied as is and all use is at your own risk. The authors disclaim all warranties of any kind, either express or implied, as to the softwares, including, but not limited to, implied warranties of fitness for a particular purpose, merchantability or non - infringement of proprietary rights. Neither this agreement nor any documentation furnished under it is intended to express or imply any warranty that the operation of the software will be error - free. Under no circumstances shall the authors of the softwares provided here be liable to any user for direct, indirect, incidental, consequential, special, or exemplary damages, arising from the software, or user' s use or misuse of the softwares. Such limitation of liability shall apply whether the damages arise from the use or misuse of the data provided or errors of the software
