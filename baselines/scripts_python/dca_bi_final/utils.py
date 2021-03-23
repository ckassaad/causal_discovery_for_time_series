import torch
import numpy as np;
from torch.autograd import Variable
from graph import *;
import time
from scipy.io import loadmat
from sklearn.preprocessing import StandardScaler

def normal_std(x):
    return x.std() * np.sqrt((len(x) - 1.)/(len(x)))

class Data_utility(object):
    # train and valid is the ratio of training set and validation set. test = 1 - train - valid
    def __init__(self, args):
        self.cuda = args.cuda
        self.model = args.model
        self.P = args.window
        self.h = args.horizon
        self.y_dim = args.y_dim        
        self.pre_win = args.pre_win 
        
        self.rawdat = loadmat(args.data)['expression']
        print('data shape', self.rawdat.shape);

        self.dat = np.zeros(self.rawdat.shape);
        self.n, self.m = self.dat.shape
        if args.normalize == 1:            
            self._normalized();

        self._bi_split(args.train, args.valid)
 

    def _normalized(self):
        
        for i in range(self.m):
            Mean = np.mean(self.rawdat[:,i])
            Std = np.std(self.rawdat[:,i])
            
            self.dat[:,i] = (self.rawdat[:,i] - Mean) / Std
        
    def _bi_split(self, train, valid):

        bi_w = 2*self.P + self.pre_win
        n = len(self.dat)
        num_biwin = n - bi_w + 1
        
        X_bi = torch.zeros((num_biwin, bi_w, self.m))
        for i in range(num_biwin):
            start = i
            end = i + bi_w
            X_bi[i, :, :] = torch.from_numpy(self.dat[start : end, :])
            
        index_bi = torch.randperm(len(X_bi))
        X_bi = X_bi[index_bi, :, :]
        
        train_set = range(0, int(train * num_biwin))
        valid_set = range(int(train * num_biwin), int((train + valid) * num_biwin))
        test_set = range(int((train + valid) * num_biwin), num_biwin)
        
        self.train_bi = X_bi[train_set, :, :]
        self.valid_bi = X_bi[valid_set, :, :]
        self.test_bi = X_bi[test_set, :, :]


    def get_batches_bi(self, data, batch_size, shuffle = False):
        
        length = len(data)
        index= torch.LongTensor(range(length))
        start_idx = 0
        
        while (start_idx < length):
            end_idx = min(length, start_idx + batch_size)
            excerpt = index[start_idx : end_idx]
            data_batch = data[excerpt]
            forward_idx  = range(0, self.P)
            Y_idx = range(self.P, self.P + self.pre_win)
            backward_idx = np.arange(data_batch.shape[1], self.P + self.pre_win, -1) - 1
            X_fw = data_batch[:, forward_idx, :] 
            Y = data_batch[:, Y_idx, :]
            X_bw = data_batch[:, backward_idx, :]
            
            if (self.cuda):
                X_fw = X_fw.cuda()
                Y = Y.cuda()
                X_bw = X_bw.cuda()
                
            model_inputs_fw = [Variable(X_fw)]
            model_inputs_bw = [Variable(X_bw)]

            data_out = [model_inputs_fw, Variable(Y), model_inputs_bw]
            yield data_out
            start_idx += batch_size            
            
            
