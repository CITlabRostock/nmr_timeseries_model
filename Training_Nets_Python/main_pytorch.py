import torch
import random as rd
import numpy as np
import torch.nn.functional as func
import matplotlib.pyplot as plt
from torch import nn
from torch import optim
from torch.utils.data import Dataset, DataLoader
from scipy.special import voigt_profile

device = (
    "cuda"
    if torch.cuda.is_available()
    else "mps"
    if torch.backends.mps.is_available()
    else "cpu"
)
print(f"Using {device} device")

# Global Variables, need to be defined and named that way
numpeaks = 1
x = np.linspace(0, 10, 1000)
function_type = 'Pseudo-Voigt'
trainlen = 15000
validlen = 1500
ModelName = 'MATLAB_net_' + str(numpeaks) + 'p'
testing = 1
epochs = 25
bsize = 25
lrate = 0.001
global_loss = 2  # If > 0, will try and account for global loss, if 2 mixed with loss in parameter space
penalty = 0.1
prob_lam = 0
prob_hwid = 0

# Defines pseudo-Voigt function
def gl(x1, params):
    y2 = np.zeros(len(x1))
    logn = -1 * np.log(2)
    cent = params[0]
    heig = params[1]
    hwid = params[2]
    if len(params) == 3:
        lamb = 0
    else:
        lamb = params[3]
    for ind in range(len(x1)):
        y2[ind] = heig * (lamb * np.exp((logn / (hwid ** 2)) * ((x1[ind] - cent) ** 2)) + (1 - lamb) * (hwid ** 2) / (
                (x1[ind] - cent) ** 2 + (hwid ** 2)))
    return y2


def voigt(x1, params):
    cent = params[0]
    heig = params[1]
    gamm = params[2]
    sigm = params[3]
    x1 = x1 - cent  # Shift whole x-axis, the Voigt-function is always centered around 0
    heig = heig / voigt_profile(0, sigm, gamm)
    y2 = voigt_profile(x1, sigm, gamm) * heig
    return y2


def make_spectrum(npeaks, ppms, plamb=0, phwid=0):
    params = np.zeros((4, npeaks))
    drow = np.zeros((1, len(ppms)))
    for ind in range(npeaks):
        params[0, ind] = rd.uniform(0, 1)  # Every parameter is in [0,1] (unscaled), center
        params[1, ind] = rd.uniform(0.1, 1)  # Height of the peak should not be too small, height
        params[2, ind] = rd.uniform(0.0001, 3)  # Half width up to 3 for a test
        params[3, ind] = rd.uniform(0.0001, 1)  # Lambda
    params = params[:, params[0].argsort()]  # Sorts whole matrix by center value (first row)
    scales = [10, 1, 0.5, 0.2]
    if function_type == 'Pseudo-Voigt':
        for ind2 in range(npeaks):
            if rd.uniform(0, 1) > 1 - phwid:
                params[2, ind2] = params[2, ind2] / 4  # Create thin peaks, only for pV-curves
            if rd.uniform(0, 1) > 1 - plamb:
                params[3, ind2] = params[3, ind2] * 5  # Create large lambdas for better recognition (hopefully)
            drow += gl(ppms,
                       [params[0, ind2] * scales[0],
                        params[1, ind2] * scales[1],
                        params[2, ind2] * scales[2],
                        params[3, ind2] * scales[3]]
                       )
    elif function_type == 'Voigt':
        for ind3 in range(npeaks):
            drow += voigt(ppms, [params[0, ind3] * 10, params[1, ind3], params[2, ind3], params[3, ind3]])
    drow += + np.random.normal(0, 0.005, len(ppms))
    params = params.flatten()
    params = params.reshape((1, -1))
    drow = drow.reshape((1, -1))
    return params, drow

# Creates spectrum from a set of parameters
def params2spec(x2, params):
    params = params.reshape(4, -1)
    y2 = np.zeros(x2.size)
    for ind4 in range(params.size(dim=1)):
        if function_type == 'Pseudo-Voigt':
            y2 += gl(x2, [params[0, ind4] * 10, params[1, ind4], params[2, ind4] / 2, params[2, ind4] / 5])
        elif function_type == 'Voigt':
            y2 += voigt(x2, [params[0, ind4] * 10, params[1, ind4], params[2, ind4], params[2, ind4]])
    return y2

 # Creates spectrum from tensor of parameters
def tensor2spec(x2, bparams):
    bparams = bparams.reshape(-1, 4, numpeaks)
    xt = torch.from_numpy(x2).float().to(device)
    yt = torch.zeros(bparams.size(dim=0), x2.size).float().to(device)
    logn = -np.log(2)
    for indb in range(bparams.size(dim=0)):
        for indp in range(bparams.size(dim=2)):
            [c, h, w, lamb] = bparams[indb, :, indp]
            c = c * 10
            w = w/2
            lamb = lamb/5
            xtrafo = ((xt - c) / w)
            xtrafo = xtrafo*xtrafo
            gtensor = torch.exp(logn * xtrafo)
            ltensor = 1 / (1 + xtrafo)
            gltensor = h * (lamb * gtensor + (1 - lamb) * ltensor)
            yt[indb] += gltensor
    return yt


class SpectraSet(Dataset):
    def __init__(self, num_samples):
        self.data = []
        perct = 5
        modul = int(np.round(num_samples / (100 / perct)))
        print('Creating Dataset of ' + str(num_samples) + ' spectra. Progress:')
        for ind5 in range(num_samples):
            paras, datarow = make_spectrum(numpeaks, x, prob_lam, prob_hwid)
            self.data.append([datarow, paras])
            if (ind5 + 1) % modul == 0:
                print(str((int(ind5 / modul) + 1) * perct) + ' %')
        print('Successfully created dataset of ' + str(num_samples) + ' spectra consisting of ' + str(numpeaks) + ' ' +
              function_type + ' functions')
        self.spec_len = len(x)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        vector, label = self.data[idx]
        vector = torch.tensor(vector, dtype=torch.float)  # Must define as float for convlayers need floats
        label = torch.tensor(label, dtype=torch.float)
        return vector, label


train_data = SpectraSet(trainlen)
valid_data = SpectraSet(validlen)
train_dataloader = DataLoader(train_data, batch_size=bsize, shuffle=True)
valid_dataloader = DataLoader(valid_data, batch_size=bsize, shuffle=True)


class myCNN(nn.Module):
    def __init__(self, output_dim):
        super().__init__()
        self.c1 = nn.Conv1d(1, 2, 32)  # Batch x hat shape [batch_size, # channels=1, len(vec)]
        self.c2 = nn.Conv1d(2, 2, 16)  # Conv1d hat args (# in_channels, #out_channels, kernel_size)
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


mymodel = myCNN(4 * numpeaks).to(device)  # Ensures GPU usage
Loss_fn = nn.MSELoss(reduction='sum')  # Doesn't average
adam_optimizer = optim.Adam(mymodel.parameters(), lr=lrate)


# Training routine
def train_loop(dataloader, model, loss_fn, optimizer):
    size = len(dataloader.dataset)
    model.train()  # Defines training mode
    if global_loss > 0 and function_type == 'Pseudo-Voigt':
        for batch, (inputs, targets) in enumerate(dataloader):
            inputs, targets = inputs.to(device), targets.to(device)

            # Compute prediction error
            pred = model(inputs)
            pred_spec = tensor2spec(x, pred)
            if global_loss > 1:
                loss1 = loss_fn(pred_spec, inputs[:, 0, :])
                loss2 = loss_fn(pred, targets)
                loss = penalty*loss1 + loss2
            else:
                loss = loss_fn(pred_spec, inputs[:, 0, :])
            loss /= bsize

            # Backpropagation
            loss.backward()
            optimizer.step()  # Optimizer gets parameters, which possess gradients
            optimizer.zero_grad()

            if (batch + 1) % 100 == 0:
                train_loss, current = loss.item(), (batch + 1) * len(inputs)
                print(f"loss: {train_loss:>7f}  [{current:>5d}/{size:>5d}]")
            if torch.isnan(loss):
                print('Loss is NaN!')
    else:
        for batch, (inputs, targets) in enumerate(dataloader):
            inputs, targets = inputs.to(device), targets.to(device)

            # Compute prediction error
            pred = model(inputs)
            loss = loss_fn(pred, targets)
            loss /= bsize
            # Backpropagation
            loss.backward()
            optimizer.step()  # Optimizer gets parameters, which possess gradients
            optimizer.zero_grad()

            if (batch + 1) % 100 == 0:
                train_loss, current = loss.item(), (batch + 1) * len(inputs)
                print(f"loss: {train_loss:>7f}  [{current:>5d}/{size:>5d}]")
    train_loss = loss.item()
    return train_loss


def valid_loop(dataloader, model, loss_fn):
    # Set the model to evaluation mode - important for batch normalization and dropout layers
    model.eval()
    num_batches = len(dataloader)
    valid_loss = 0

    # Evaluating the model with torch.no_grad() ensures that no gradients are computed during test mode
    # also serves to reduce unnecessary gradient computations and memory usage for tensors with requires_grad=True
    if global_loss > 0 and function_type == 'Pseudo-Voigt':
        with torch.no_grad():
            for inputs, targets in dataloader:
                pred = model(inputs)
                pred_spec = tensor2spec(x, pred)
            if global_loss > 1:
                loss1 = loss_fn(pred_spec, inputs[:, 0, :])
                loss2 = loss_fn(pred, targets)
                loss = penalty*loss1 + loss2
            else:
                loss = loss_fn(pred_spec, inputs[:, 0, :])
            loss /= bsize
            valid_loss += loss
        print(f"Avg test loss: {valid_loss:>8f} \n")
        return valid_loss
    else:
        with torch.no_grad():
            for inputs, targets in dataloader:
                pred = model(inputs)
                valid_loss += loss_fn(pred, targets).item()

        valid_loss /= (num_batches*bsize)
        print(f"Avg test loss: {valid_loss:>8f} \n")
        return valid_loss


plotloss = np.zeros((2, epochs))
for j in range(epochs):
    print(f"Epoch {j + 1}\n-------------------------------")
    t_loss = train_loop(train_dataloader, mymodel, Loss_fn, adam_optimizer)
    v_loss = valid_loop(valid_dataloader, mymodel, Loss_fn)
    plotloss[0][j] = t_loss
    plotloss[1][j] = v_loss
print("Done!")

# Plotting results

plt.plot(range(epochs), plotloss[0], range(epochs), plotloss[1])
plt.ylabel('Losses')
plt.xlabel('Epoch')
plt.legend(['Train_Loss', 'Val_loss'], loc='upper right')
plt.savefig('Bilder/pytorch/' + ModelName + '.png')
plt.show()
plt.close()

torch.save(mymodel.state_dict(), 'models/pytorch/' + ModelName + '_SD.pt', _use_new_zipfile_serialization=False) # Saving the networks

if testing == 1:
    with torch.no_grad():
        inp, tar = next(iter(valid_dataloader))
        inp = inp.to(device)
        tar = tar.to(device)
        result = mymodel(inp)
        for i in range(min(10, inp.size(dim=0))):
            y = params2spec(x, result[i, 0])
            plt.plot(x, y, x, inp[i, 0])
            plt.ylabel('Amplitude')
            plt.xlabel('Chem. shift')
            plt.legend(['Prediction', 'Inputs'], loc='upper right')
            plt.show()
            plt.close()
            
