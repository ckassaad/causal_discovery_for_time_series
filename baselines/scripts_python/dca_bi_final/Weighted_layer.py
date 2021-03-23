import torch
import torch.nn as nn
import pdb
import torch.nn.functional as F
from torch.nn.parameter import Parameter
import numpy as np
from sklearn.preprocessing import normalize
## granger full rank

class WL_Model(nn.Module):
    def __init__(self):
        super(WL_Model, self).__init__()
        self.weight = nn.Parameter(torch.ones([1]))
        #pdb.set_trace()
        #self.weight.data = torch.nn.Parameter(torch.from_numpy(data.GTs).float()).reshape(lowrank,1);
        
        
    def forward(self, x):
        #l = x.shape[1]
        y = torch.mul(x, self.weight)
        return y
        


    

