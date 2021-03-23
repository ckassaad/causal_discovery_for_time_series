import torch
import torch.nn as nn
import pdb
import torch.nn.functional as F
from torch.nn.parameter import Parameter
import FinalRowWise_layer, SingularValue_layer, BiDirection_layer
import numpy as np
from sklearn.preprocessing import normalize
import matplotlib.pyplot as plt
from torch.autograd import Variable
## granger full rank
class Model(nn.Module):
    def __init__(self, args, data):
        super(Model, self).__init__()
        self.use_cuda = args.cuda

        self.m = data.m
        self.w = args.window

        
        self.batch_size = args.batch_size
        self.lowrank = args.lowrank

        self.pre_win = args.pre_win 
        self.p_list = args.p_list
        self.y_dim = args.y_dim
        self.p_allsum = np.sum(self.p_list)
        self.len_p_list = len(self.p_list)
        self.compress_p_list = args.compress_p_list
        self.len_compress_p_list = len(self.compress_p_list)

        self.sparse_label = []
        self.orthgonal_label = []
        
        self.linears = [(nn.Linear(self.w, self.p_list[0]))] #w->hid
        self.sparse_label.append(0)
        self.orthgonal_label.append(1)

        if self.len_p_list>1:
            for p_i in np.arange(1, self.len_p_list):
                self.linears.append((nn.Linear(self.p_list[p_i-1], self.p_list[p_i]))) #w->hid
                self.sparse_label.append(0)
                self.orthgonal_label.append(1)

        
        ## graph layers
        for p_i in np.arange(0, self.len_p_list):
            self.linears.append(BiDirection_layer.BD_Model(self.m, self.lowrank))
            self.sparse_label.append(1)
            self.orthgonal_label.append(1)
            
            self.linears.append(nn.BatchNorm1d(self.p_list[-1]))
            self.sparse_label.append(0)
            self.orthgonal_label.append(0)
            
            self.linears.append(nn.BatchNorm1d(self.p_list[-1]))
            self.sparse_label.append(0)
            self.orthgonal_label.append(0)
            
            
            
            
            self.linears.append(SingularValue_layer.SV_Model(data, self.lowrank, self.use_cuda))
            self.sparse_label.append(0)
            self.orthgonal_label.append(0)
            
            self.linears.append(nn.BatchNorm1d(self.p_list[-1]))
            self.sparse_label.append(0)
            self.orthgonal_label.append(0)
            
            self.linears.append(nn.BatchNorm1d(self.p_list[-1]))
            self.sparse_label.append(0)
            self.orthgonal_label.append(0)
            
    
            self.linears.append(BiDirection_layer.BD_Model(self.lowrank, self.m))
            self.sparse_label.append(1)
            self.orthgonal_label.append(1)
            
            self.linears.append(nn.BatchNorm1d(self.p_list[-1]))
            self.sparse_label.append(0)
            self.orthgonal_label.append(0)
            
            self.linears.append(nn.BatchNorm1d(self.p_list[-1]))
            self.sparse_label.append(0)
            self.orthgonal_label.append(0)
            
                    
        if self.len_compress_p_list>0:
            self.linears.append(nn.Linear(self.p_allsum, self.compress_p_list[0]))
            self.sparse_label.append(0)
            self.orthgonal_label.append(1)
            for p_j in np.arange(1, self.len_compress_p_list):
                self.linears.append(nn.Linear(self.compress_p_list[p_j-1], self.compress_p_list[p_j]))
                self.sparse_label.append(0)
                self.orthgonal_label.append(1)
                
        if self.len_compress_p_list>0:
            self.linears.append(nn.Linear(self.p_allsum, self.compress_p_list[0]))
            self.sparse_label.append(0)
            self.orthgonal_label.append(1)
            for p_j in np.arange(1, self.len_compress_p_list):
                self.linears.append(nn.Linear(self.compress_p_list[p_j-1], self.compress_p_list[p_j]))
                self.sparse_label.append(0)
                self.orthgonal_label.append(1)                
                
                
                
          
        self.linears.append(FinalRowWise_layer.FR_Model(args, 2*self.compress_p_list[-1])) #k->k  
        self.sparse_label.append(1)
        self.orthgonal_label.append(0)
        
        self.linears = nn.ModuleList(self.linears)
        self.dropout = nn.Dropout(args.dropout)


    def forward(self, x_fw, x_bw):
        x_fw = x_fw[0]
        x_bw = x_bw[0]

        x_fw = x_fw.transpose(2, 1).contiguous()
        x_fw = self.dropout(x_fw)
        x_fw_org = x_fw
        x_fw_p = []
                
        if self.p_list[0] > self.w:
            padding = nn.ConstantPad2d((0, self.p_list[0]-self.w, 0, 0), 0)
            x_0n = padding(x_fw_org)       
        x_0 = x_fw_org
        for layer_i in range(self.len_p_list):  
            x_i = self.linears[layer_i + 0](x_0)
            x_i = F.relu(x_i + x_0n)
            x_0n = x_i
            x_0 = x_i
            x_i = self.dropout(x_i)
            x_fw_p.append(x_i)
            
            
        x_bw = x_bw.transpose(2, 1).contiguous()
        x_bw = self.dropout(x_bw)
        x_bw_org = x_bw
        x_bw_p = []
                
        if self.p_list[0] > self.w:
            padding = nn.ConstantPad2d((0, self.p_list[0]-self.w, 0, 0), 0)
            x_0n = padding(x_bw_org)       
        x_0 = x_bw_org
        for layer_i in range(self.len_p_list):  
            x_i = self.linears[layer_i + 0](x_0)
            x_i = F.relu(x_i + x_0n)
            x_0n = x_i
            x_0 = x_i
            x_i = self.dropout(x_i)
            x_bw_p.append(x_i)            
            
            
        x_p_m_fw = []  
        for layer_i in range(self.len_p_list):
            
            x_sp =  x_fw_p[layer_i].transpose(2,1).contiguous()
            
            x_sp = self.linears[self.len_p_list + layer_i*9 + 0](x_sp, direction = 'forward')
            x_sp = self.linears[self.len_p_list + layer_i*9 + 1](x_sp)
            x_sp = F.tanh(x_sp/5.)
            x_sp = self.dropout(x_sp)
            
            x_sp = self.linears[self.len_p_list + layer_i*9 + 3](x_sp)
            x_sp = self.linears[self.len_p_list + layer_i*9 + 4](x_sp)
            x_sp = self.dropout(x_sp)
            
            x_sp = self.linears[self.len_p_list + layer_i*9 + 6](x_sp, direction = 'forward')
            x_sp = self.linears[self.len_p_list + layer_i*9 + 7](x_sp)
            x_sp = F.selu(x_sp/1.)
            x_sp = self.dropout(x_sp)
            x_sp = x_sp.transpose(2, 1).contiguous()
            x_p_m_fw.append(x_sp)            
            
            
        x_p_m_bw = []  
        for layer_i in range(self.len_p_list):
            x_sp =  x_bw_p[layer_i].transpose(2,1).contiguous()
            
            x_sp = self.linears[self.len_p_list + layer_i*9 + 6](x_sp, direction = 'backward')
            x_sp = self.linears[self.len_p_list + layer_i*9 + 8](x_sp)
            x_sp = F.tanh(x_sp/5.)
            x_sp = self.dropout(x_sp)
            
            x_sp = self.linears[self.len_p_list + layer_i*9 + 3](x_sp)
            x_sp = self.linears[self.len_p_list + layer_i*9 + 5](x_sp)
            x_sp = self.dropout(x_sp)
            
            x_sp = self.linears[self.len_p_list + layer_i*9 + 0](x_sp, direction = 'backward')
            x_sp = self.linears[self.len_p_list + layer_i*9 + 2](x_sp)
            x_sp = F.selu(x_sp/1.)
            x_sp = self.dropout(x_sp)
            x_sp = x_sp.transpose(2, 1).contiguous()
            x_p_m_bw.append(x_sp)              
            
            
        x_p_m_fw = torch.cat(x_p_m_fw, dim = 2)
        x_p_m_fw = x_p_m_fw[:, 0:self.y_dim, :]
            
            
        if self.len_compress_p_list > 0:
            for p_j in range(self.len_compress_p_list): 
                x_p_m_fw = self.linears[self.len_p_list + self.len_p_list * 9 + p_j](x_p_m_fw) 
                x_p_m_fw = F.tanh(x_p_m_fw/5.)
                x_p_m_fw = self.dropout(x_p_m_fw)
                                
        x_p_m_bw = torch.cat(x_p_m_bw, dim = 2)
        x_p_m_bw = x_p_m_bw[:, 0:self.y_dim, :]
            
            
        if self.len_compress_p_list > 0:
            for p_j in range(self.len_compress_p_list): 
                x_p_m_bw = self.linears[self.len_p_list + self.len_p_list * 9 + self.len_compress_p_list + p_j](x_p_m_bw) 
                x_p_m_bw = F.tanh(x_p_m_bw/5.)
                x_p_m_bw = self.dropout(x_p_m_bw)
         
        x_p_m = torch.cat((x_p_m_fw, x_p_m_bw), 2)
        final_y = self.linears[-1](x_p_m)

        return final_y      
 


    
    def predict_relationship(self):
        
        #weight_help1, weight_help2 = self.linears[-1].get_pi_weight()
        
        CGraph_list1 = []
        CGraph_list2 = []
        G_1 = np.zeros((self.m,self.m))
        G_2 = np.zeros((self.m,self.m))
        G_3 = np.zeros((self.m,self.m))
        G_4 = np.zeros((self.m,self.m))
        
                
        for layer_i in range(self.len_p_list):
            
            A = self.linears[self.len_p_list + layer_i * 9 + 0].weight.transpose(0,1).cpu().detach().numpy()
            B = np.diag(self.linears[self.len_p_list + layer_i * 9 + 3].weight.transpose(0,1).detach().cpu().numpy().ravel())
            C = self.linears[self.len_p_list + layer_i*9 + 6].weight.transpose(0,1).cpu().detach().numpy()

            CGraph1 = np.abs(np.dot(np.dot(A, B), C))
            CGraph1[range(self.m), range(self.m)] = 0    
            CGraph_list1.append(CGraph1)
            
            A = np.abs(A) 
            B = np.abs(B) 
            C = np.abs(C) 

            CGraph2 = np.abs(np.dot(np.dot(A, B), C))
            CGraph2[range(self.m), range(self.m)] = 0    
            CGraph_list2.append(CGraph2)   

            
            G_1 = np.add(G_1, CGraph1)
            G_2 = np.add(G_2, CGraph2)
            
            G_3 = G_1
            G_4 = G_2
            
        G_1[range(self.m), range(self.m)] = 0 
        G_2[range(self.m), range(self.m)] = 0 
        G_3[range(self.m), range(self.m)] = 0 
        G_4[range(self.m), range(self.m)] = 0 
        
           
        return G_1, G_3, G_2, G_4