#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 12 14:52:36 2019

@author: chenxiaoxu
"""

import torch
import torch.nn as nn
import pdb
import torch.nn.functional as F
from torch.nn import init
from torch.nn.parameter import Parameter
import numpy as np
import math

class BD_Model(nn.Module):
    def __init__(self, in_features, out_features):
        super(BD_Model, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.weight = Parameter(torch.Tensor(out_features, in_features))
        init.kaiming_uniform_(self.weight, a = math.sqrt(5))
        
        
    def forward(self, inputs, direction = 'forward'):
        if direction == 'forward':
            return F.linear(inputs, self.weight)
        elif direction == 'backward':
            return F.linear(inputs, self.weight.t())
        
        
        