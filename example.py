import functions as deepshm 
# when you are actually using the package you will have just 
# import deepshm
import numpy as np

###########################
# generate data
###########################
# model's input must have shape (n, 4, k, 1), 
# where n is the number of training/test sequences and k is the size of k-mer
# model's target must have shape (n, 4) for 'substitution' and 'weighted_sub' tasks or (n) for 'mut_freq' task
n_train = 1000
n_test = 100
k = 15 # size of the kmer 
# if you want to use one of the provided models k has to be 5, 9, 15 or 21
# if you want to train your own model k can be any positive integer

train_inputs = np.zeros((n_train, 4, k, 1))             # create an array for one-hot encoded train inputs
train_seqs = np.random.randint(4, size=(n_train, k))    # generate n_train sequences of length k. 0 = A, 1 = C, 2 = G, 3 = T
# transform to the inpur format
for i in range(n_train):
    for j in range(k):
        train_inputs[i,train_seqs[i,j],j,0] = 1
train_targets = np.random.rand(n_train, 4)              # should have shape (n_train, 4) for 'substitution' and 'weighted_sub' tasks
                                                        # or (n_train) for 'mut_freq' task
test_inputs = np.zeros((n_test, 4, k, 1))               # create an array for one-hot encoded test inputs
test_seqs = np.random.randint(4, size=(n_test, k))      # generate n_test sequences of length k. 0 = A, 1 = C, 2 = G, 3 = T
# transform to the inpur format
for i in range(n_test):
    for j in range(k):
        test_inputs[i,test_seqs[i,j],j,0] = 1
test_targets = np.random.rand(n_test, 4)            # should have shape (n_test, 4) for 'substitution' and 'weighted_sub' tasks
                                                    # or (n_test) for 'mut_freq' task
# combine train and test data in one list
data = [train_inputs, train_targets, test_inputs, test_targets]

###########################
# initialize model parameters
###########################
#set up network
channels_conv = [256]   # the length of the list determines the number of convolutional layers
kernel_sizes = [[4,9]]  # must be the same length with channels_conv
pooling = [1]           # must be the same length with channels_conv
dropout_conv = [0.4]    # must be the same length with channels_conv
channels_fc = [256,8]   # the length of the list determines the number of dense layers
dropout_fc = [0.1,0.5]  # must be the same length with channels_fc
#set up training parameters
n_epochs = 150
learning_rate = np.power(0.5, 14)
mini_batch_size =  16
task = 'weighted_sub'
#path for saving final model
save_path = 'current_model_test.h5'

###########################
# train model
###########################
training_res = deepshm.train_model(
    data,
    channels_conv,
    kernel_sizes,
    pooling,
    dropout_conv,
    channels_fc,
    dropout_fc,
    task,
    save_path=save_path,
    n_epochs=n_epochs,
    learning_rate=learning_rate,
    mini_batch_size=mini_batch_size,
)
print('Model results')
print(training_res)