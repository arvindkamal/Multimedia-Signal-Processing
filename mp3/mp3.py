#!/usr/bin/env python
# coding: utf-8

# # Machine Problem 3: NumPy CNN

# ### Import Libraries

# In[1]:


import numpy as np
import os
import time
from scipy import signal
from imageio import imread
from random import shuffle
from matplotlib import pyplot as plt

import random
get_ipython().magic(u'matplotlib inline')


# ### Preprocessing Functions

# In[2]:


# load_images
    # Read in images and makes a list for each set in the form: [images, labels]
    # images: np array with dims [N x img_height x img width x num_channels]
    # labels: np array with dims [N x 1]. elephant = 0, lionfish = 1
    #
    # Returns:  train_set: The list [train_images, train_labels]
    #           val_set: The list [val_images, val_labels] 

def load_images():
    
    sets = ['train', 'val']
    
    data_sets = []
    for dset in sets:
        img_path = './bin_dataset/' + dset + '/ele'
        ele_list = [imread(os.path.join(img_path, img))/255 for img in os.listdir(img_path)]

        img_path = './bin_dataset/' + dset + '/lio'
        lio_list = [imread(os.path.join(img_path, img))/255 for img in os.listdir(img_path)]

        set_images = np.stack(ele_list + lio_list)
        N = set_images.shape[0]
        labels = np.ones((N,1))
        labels[0:int(N/2)] = 0
        data_sets.append([set_images, labels])

    train_set, val_set = data_sets

    print("Loaded", len(train_set[0]), "training images")
    print("Loaded", len(val_set[0]), "validation images")
    
    return train_set, val_set


# batchify
    # Inputs:    train_set: List containing images and labels
    #            batch size: The desired size of each batch
    #
    # Returns:   image_batches: A list of shuffled training image batches, each with size batch_size
    #            label_batches: A list of shuffled training label batches, each with size batch_size 

def batchify(train_set, batch_size):
    
    # YOUR CODE HERE
    l = len(train_set[0])
    idx = np.arange(l)
    shuffle(idx)
    
    train_data = train_set[0][idx]
    train_label = train_set[1][idx]
    
    image_batches = []
    label_batches = []
    
    for i in range(0, l, batch_size):
        image_batches.append(train_data[i:min(i+batch_size, l)])
        label_batches.append(train_label[1][i:min(i+batch_size, l)])
    
    return image_batches, label_batches


# ### Network Functions

# #### Activation Functions

# In[15]:


# relu
    # Inputs:   x: Multi-dimensional array with size N along the first axis
    # 
    # Returns:  out: Multi-dimensional array with same size of x 

def relu(x):
    
    # YOUR CODE HERE
    out = np.maximum(x, 0)
    
    return out


# sigmoid
    # Inputs:    x: Multi-dimensional array with size N along the first axis
    # 
    # Returns:   out: Multi-dimensional array with same size of x 

def sigmoid(x):
    
    # YOUR CODE HERE
    out = 1/(1+np.exp(-x))
    
    return out


# unit_step
    # Inputs:    x: Multi-dimensional array with size N along the first axis 
    # 
    # Returns:   out: Multi-dimensional array with same size of x 

def unit_step(x):
    
    # YOUR CODE HERE
    out = np.heavyside(x, 1)
    
    return out 


# #### Layer Functions

# In[12]:


# convolve2D
    # Inputs:    X: [N x height x width x num_channels]
    #            filters: [num_filters x filter_height x filter_width x num_input_channels]
    # 
    # Returns:   Xc: output array by convoling X and filters. [N x output_height x output_width x num_filters]

def convolve2D(X0, filters):
   
    N, X0_len, _, num_ch = X0.shape
    num_out_ch, filter_len, _, _ = filters.shape
    F0_side = X0_len - filter_len + 1
    
    F0 = np.zeros((N, F0_side, F0_side, num_out_ch))
    
    for n in range(N):
        for o_ch in range(num_out_ch):
            for ch in range(num_ch):
                # YOUR CODE HERE
                
                for f0_row in range(F0_side):
                    for f0_col in range(F0_side):
#                         conv_out_1_pixel = 0
#                         for fl_row in range(filter_len):
#                             for fl_col in range(filter_len):
#                                 conv_out_1_pixel += X0[n, f0_row, f0_col, ch]*filters[o_ch, :, :, ch]
                        F0[n, f0_row, f0_col, o_ch] = (X0[n, f0_row:(f0_row+filter_len), f0_col:(f0_col + filter_len), ch]*filters[o_ch, :, :, ch]).sum()
    return F0


# maxPool
    # Inputs:    R0: [N x height x width x num_channels]
    #            mp_len: size of max pool window, also the stride for this MP
    # 
    # Returns:   p_out: output of pooling R0. [N x output_height x output_width x num_channels]
    #            R0_mask: A binary mask with the same size as R0. Indicates which index was chosen to be the max
    #            for each max pool window. This will be used for backpropagation.

def maxPool(R0, mp_len):

    N, R0_len, _, num_ch = R0.shape
    p_out_len = int((R0_len-mp_len)/mp_len + 1)

    R0_mask = np.zeros(R0.shape)
    p_out = np.zeros((N, p_out_len, p_out_len, num_ch))
    
    for n in range(N):
        for ch in range(num_ch):
            for row in range(p_out_len): 
                for col in range(p_out_len):
                    # YOUR CODE HERE
                    
                    pool = R0[row*mp_len: (row+1)*mp_len, col*mp_len: (col+1)*mp_len]
                    p_out[n, row, col, ch] = pool.max()
                    
                    for R0_row in range(row*mp_len, (row+1)*mp_len):
                        for R0_col in range(col*mp_len, (col+1)*mp_len):
                            x, y = pool.argmax()/mp_len, pool.argmax()%mp_len
                            if(R0_row==x and R0_col==y):
                                R0_mask[n, R0_row, R0_col, ch] = 1
                            else:
                                R0_mask[n, R0_row, R0_col, ch] = 0

    return p_out, R0_mask


# fc
    # Inputs:    X: [N x num_input_features]
    #            W: [num_input_features x num_fc_nodes]
    # 
    # Returns:   out: Linear combination of X and W. [N x num_fc_nodes]

def fc(X, W):
    
    # YOUR CODE HERE
    out = X @ W
    return out


# #### CNN Functions

# In[13]:


# cnn_fwd
    # Inputs:    X0: batch of images. [N x img_height x img_width x num_channels]
    #            W0, W1, W2: Parameters of the CNN
    #            mp_len: the length of one side of the max pool window
    # 
    # Returns:   sig: vector containing the output for each sample. [N x 1]
    #            cache: a dict containing the relevant output layer calculations that will be
    #            used in backpropagation
    
def cnn_fwd(X0, W0, W1, W2, mp_len):
    # F0 
    # YOUR CODE HERE
    F0 = convolve2D(X0, W0)
    R0 = relu(F0)
    
    
    # X1p 
    # YOUR CODE HERE
    X1p, R0_mask = maxPool(R0, mp_len)
    
    # X1 (flatten)
    # YOUR CODE HERE
    X1 = np.zeros((X1p.shape[0], X1p.shape[1]*X1p.shape[2]))
    for ch in range(X1p.shape[3]):
        X1 += X1p[N, :, :, ch].flatten()
    
    
    # FC Layers
    # YOUR CODE HERE
    F1 = fc(X1, W1)
    X2 = relu(F1)
    F2 = fc(X2, W2)
    
    # Output
    # YOUR CODE HERE
    sig = sigmoid(F2)
    
    # Save outputs of functions for backward pass
    cache = {
        "F0":F0,
        "R0":R0,
        "X1p":X1p,
        "R0m":R0_mask,
        "X1":X1,
        "F1":F1,
        "X2":X2,
        "F2":F2      
    }
    
    return sig, cache


# loss
    # Inputs:    sig: vector containing the CNN output for each sample. [N x 1]
    #            Y: vector containing the ground truth label for each sample. [N x 1]
    # 
    # Returns:   L: Loss/error criterion for the model. 

def loss(sig, Y):
    # YOUR CODE HERE
    L = -np.sum(np.dot(Y, np.log(sig)))/Y.shape[0]
    return L


# ### Backprop Functions

# In[6]:


# convolve2DBwd
    # Inputs:    X0: batch of images. [N x height x width x num_channels]
    #            dL_dF0: Gradient at the output of the conv layer. 
    # 
    # Returns:   dL_dW0. gradient of loss L wrt W0. Same size as W0

def convolve2DBwd(X0, dL_dF0):
    
    N, X0_len, _, num_ch = X0.shape
    _, dL_dF0_len, _, num_out_ch  = dL_dF0.shape
    filter_len = X0_len - dL_dF0_len + 1
    
    dL_dW0 = np.zeros((num_out_ch, filter_len, filter_len, num_ch))
    
    for n in range(N):
        for o_ch in range(num_out_ch):
            for ch in range(num_ch):
                # YOUR CODE HERE 
    
    return dL_dW0


# maxPoolBwd
    # Inputs:    dL_dX1p: Gradient at the output of the MaxPool layer
    #            R0_mask: A binary mask with the same size as R0. Defined in maxPool
    #            mp_len: the length of one side of the max pool window
    # 
    # Returns:   dL_dR0: Gradient at the output of ReLu
    
def maxPoolBwd(dL_dX1p, R0_mask,  mp_len):
    
    N, H, W, C = R0_mask.shape
    N, dH, dW, C = dL_dX1p.shape
    
    dL_dR0 = np.zeros(R0_mask.shape)
    
    for n in range(N):
        for ch in range(C):
            for row in range(dH):
                for col in range(dW):
                    # YOUR CODE HERE
                    
    return dL_dR0


# dL_dW2
    # Inputs:    Y: vector containing the ground truth label for each sample. [N x 1]
    #            cache: a dict containing the relevant output layer calculations 
    # 
    # Returns:   dL_dW2: Gradient of the Loss wrt W2
    
def dL_dW2(Y, cache):
   
    # YOUR CODE HERE
    
    return dL_dW2


# dL_dW1
    # Inputs:    Y: vector containing the ground truth label for each sample. [N x 1]
    #            W2: Weight matrix for the second FC layer
    #            cache: a dict containing the relevant output layer calculations 
    # 
    # Returns:   dL_dW1: Gradient of the Loss wrt W1
    
def dL_dW1(Y, W2, cache):
    
    # YOUR CODE HERE
    
    return dL_dW1


# dL_dW0
    # Inputs:    X0: batch of images. [N x height x width x num_channels]
    #            Y: vector containing the ground truth label for each sample. [N x 1]
    #            W1: Weight matrix for the first FC layer
    #            W2: Weight matrix for the second FC layer
    #            mp_len: the length of one side of the max pool window
    #            cache: a dict containing the relevant output layer calculations 
    # 
    # Returns:   dL_dW0: Gradient of the Loss wrt W0

def dL_dW0(X0, Y, W1, W2, mp_len, cache):
    
    N, X1p_len, _, no_out_ch  = cache['X1p'].shape
    F2 = cache['F2']
    F1 = cache['F1']
    R0m = cache['R0m']
    F0 = cache['F0']
    
    #dL_dF2
    # YOUR CODE HERE
    
    #dL_dF1
    # YOUR CODE HERE
    
    #dL_dX1
    # YOUR CODE HERE
    
    # dL_dX1p (unflatten)
    # YOUR CODE HERE
    
    # dL_dR0 (unpool)
    # YOUR CODE HERE
    
    # dL_dF0 (relu_bwd)
    # YOUR CODE HERE
    
    # dL_dW0
    # YOUR CODE HERE
    
    return dL_dW0
    
        


# ### Training

# #### Load Images

# In[7]:


# Load images and scale them
# YOUR CODE HERE
train_set, val_set = load_images()


# In[314]:


print(train_set[0][0][0][0])
print(len(val_set[0]))


# In[10]:


img_batches, label_batches = batchify(train_set, batch_size) 


# In[319]:


print(len(img_batches[0][0][0]))


# #### Config

# In[9]:


# Hyperparameters
epochs = 20
lr = 0.1
batch_size = 16
filter_len = 5
num_out_ch = 3
mp_len = 12
fc_nodes = 2

# Declare weights
# YOUR CODE HERE


# In[ ]:


for i in range(epochs):
    
    # make set of batches
    # YOUR CODE HERE
    img_batches, label_batches = batchify(train_set, batch_size) 
    for b_idx in range(num_batches):
        X = img_batches[b_idx]
        Y = label_batches[b_idx]
        
        # Forward pass
        # YOUR CODE HERE
        
        # Calculate gradients
        # YOUR CODE HERE
        
        # Update gradients
        # YOUR CODE HERE
      


# ### Test Correctness of Forward and Backward Pass

# #### Forward

# In[16]:


weights = np.load('weights.npz')
W0 = weights['W0']
W1 = weights['W1']
W2 = weights['W2']

sig, _ = cnn_fwd(val_set[0], W0, W1, W2, mp_len)
train_acc = len(np.where(np.round(sig) == val_set[1])[0])/len(val_set[1])

print("train_loss:", loss(sig, val_set[1]), "train_acc:", train_acc)


# #### Backward

# In[ ]:


# Make backprop testing batch
X_bp = np.vstack([train_set[0][0:8,:,:,:], train_set[0][-9:-1,:,:,:]])
Y_bp = np.vstack([train_set[1][0:8], train_set[1][-9:-1]])

# Initialize weights to all ones
# YOUR CODE HERE

# Update weights once
# YOUR CODE HERE


print("W2 value:", np.sum(W2))
print("W1 value:", np.sum(W1))
print("W0 value:", np.sum(W0))

