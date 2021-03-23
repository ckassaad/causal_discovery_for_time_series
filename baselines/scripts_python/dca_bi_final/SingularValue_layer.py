import torch
import torch.nn as nn
import pdb
import torch.nn.functional as F
from torch.nn.parameter import Parameter
import numpy as np
from sklearn.preprocessing import normalize
## granger full rank

class SV_Model(nn.Module):
    def __init__(self, data, lowrank, use_cuda):
        super(SV_Model, self).__init__()
        self.weight = nn.Parameter(torch.ones([lowrank, 1]))

        self.use_cuda = use_cuda

    def forward(self, x):
        k = x.shape[2]        
        y = torch.Tensor(x.shape)
        
        for j in range(k):            
            tmp_new = torch.mul(x[:, :, j], self.weight[j, 0])
            y[:, :, j] = tmp_new
            
        if self.use_cuda:
            y = y.cuda()
            
        return y