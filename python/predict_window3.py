# variable "longwin" will be given by matlab
import torch
import numpy as np
from torch import nn
import torch.nn.functional as func

longwin = np.asarray(longwin)
numsplits = round(len(longwin)/1000) # The long window consists of numsplits many actual subwindows


# Definition of the CNN architecture
class myCNN(nn.Module):
    def __init__(self, output_dim):
        super().__init__()
        self.c1 = nn.Conv1d(1, 2, 32)  # Batch x has shape [batch_size, # channels=1, len(vec)]
        self.c2 = nn.Conv1d(2, 2, 16)  # 'Conv1d' has args (# in_channels, #out_channels, kernel_size)
        self.pool = nn.MaxPool1d(2, stride=2)
        self.linear_sigmoid_stack = nn.Sequential(
            nn.LazyLinear(320),
            nn.Sigmoid(),
            nn.LazyLinear(160),
            nn.Sigmoid(),
            nn.LazyLinear(80),
            nn.Sigmoid(),
            nn.LazyLinear(60),
            nn.Sigmoid(),
        )
        self.outputlayer = nn.LazyLinear(output_dim)

    def forward(self, t):
        t = self.pool(func.sigmoid(self.c1(t)))
        t = self.pool(func.sigmoid(self.c2(t)))
        t = self.linear_sigmoid_stack(t)
        t = t.reshape(t.size(dim=0), 1, -1)
        t = self.outputlayer(t)
        return t


with torch.no_grad(): # Loads CNNs from state_dicts, then predicts parameters, stores everything and gives back to matlab.
    out_d = {} # Output given back to Matlab
    for numpeaks in range(5):
        ModelName = 'MATLAB_net_' + str(numpeaks+1) + 'p'
        path = 'python/pytorch/' + ModelName + '_SD.pt'
        model = myCNN(4*(numpeaks+1))
        model.load_state_dict(torch.load(path))  # Call 'torch.load(path)' to get the dict object for 'load_state_dict'
        model.eval()  # Only load model once
        for s in range(numsplits):
            win = longwin[range(1000*s, 1000*(s+1))] # Defining current window
            win = torch.from_numpy(win)
            win = win.to(dtype=torch.float32)
            res = model(win.reshape(1, 1, -1))
            res = res.flatten()
            out_d['s' + str(s+1) + 'p' + str(numpeaks+1)] = res.numpy()