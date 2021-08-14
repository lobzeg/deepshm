import numpy as np
import h5py
import tensorflow as tf
from tensorflow import keras
from scipy import stats


SEED = 42

def one_hot(seq):
    one_hot_seq = np.zeros((4, len(seq)))
    for i in range(len(seq)):
            if seq[i] == "A":
                one_hot_seq[0][i] = 1
            if seq[i] == "C":
                one_hot_seq[1][i] = 1
            if seq[i] == "G":
                one_hot_seq[2][i] = 1
            if seq[i] == "T":
                one_hot_seq[3][i] = 1
    return one_hot_seq
 
def give_midone(test_inputs):   # provides a vector of the middle nucleotides for an array of one-hot encoded k-mers
    midone = np.zeros(test_inputs.shape[0])
    for i in range(test_inputs.shape[0]):
        for j in range(4):
            if test_inputs[i][j][int(test_inputs.shape[2]/2)][0] == 1:
                midone[i] = j
    return midone

def subs_for_corr(dset, midone):    # prepprocesses substitution predictions for correlation calculation 
                                    # deletes substitution rate corresponding to no-mutation event and makes other 3 substitution rates to add up to one
    m_f = np.zeros((dset.shape[0],1))
    new_dset = np.zeros((dset.shape[0],dset.shape[1]-1))
    for i in range(dset.shape[0]):
        h = 0
        midone_in = midone[i]
        for j in range(dset.shape[1]):
            if midone_in != j:
                new_dset[i][h] = dset[i][j]
                h+=1
        m_f[i][0] = np.sum(new_dset[i])
    subs = new_dset/m_f
    return m_f, subs

def subs_for_pred(dset, midone):    # normalizes substitution predictions 
                                    # sets up substitution rate corresponding to no-mutation event to zero and makes other 3 substitution rates to add up to one
    m_f = np.zeros((dset.shape[0],1))
    new_dset = np.zeros((dset.shape[0],dset.shape[1]))
    for i in range(dset.shape[0]):
        h = 0
        midone_in = midone[i]
        for j in range(dset.shape[1]):
            if midone_in != j:
                new_dset[i][j] = dset[i][j]
                #h+=1
        m_f[i][0] = np.sum(new_dset[i])
    subs = new_dset/m_f
    return m_f, subs

def give_corr(predictions,targets,midone,task = 'weighted_sub'): # calculates correlation between predicted values and ground truth
    if task == 'mut_freq':
        mf_corr, _ = stats.mstats.pearsonr(predictions,targets)
        subs_corr = 0
    if task == 'substitution':
        _, subs_targets = subs_for_corr(targets,midone)
        _, subs_preds = subs_for_corr(predictions,midone)
        mf_corr = 0
        subs_corr, _ = stats.mstats.pearsonr(subs_preds,subs_targets)
    if task == 'weighted_sub':
        mf_targets, subs_targets = subs_for_corr(targets,midone)
        mf_preds, subs_preds = subs_for_corr(predictions,midone)
        mf_corr, _ = stats.mstats.pearsonr(mf_preds,mf_targets)
        subs_corr, _ = stats.mstats.pearsonr(subs_preds,subs_targets)
    return mf_corr, subs_corr

def generate_model( # generates the CNN model according to parameters provided
    input_shape,    # [4*a, k, 1]
    channels_conv,  # number of channels for convolutional layers (including transition layer), length N
    kernel_sizes,   # kernel sizes for convolutional layers, length N-1
    pooling,        # kernel sizes for pooling layers following each convolutional layer (pooling[i]=0 means no pooling layer there), length N-1
    dropout_conv,   # dropout rates for convolutional layers (including transition layer), length N
    channels_fc,    # number of channels for dense layers, length M
    dropout_fc,     # dropout rates for dense layers, length M-1
    task,           # can be 'mut_freq', 'substitution', 'weighted_sub' or 'my_own'
    emb,
    padding='same',  # use 'same' padding if you want to use strand separation or 'valid' othervise
    initializer = 'he_normal',
    activation = tf.nn.elu,
):
    tf.random.set_seed(SEED) ## For reproducibility
    model = keras.Sequential()
    num_conv = len(channels_conv)
    num_fc = len(channels_fc)
    width_int = input_shape[1]
    model.add(keras.layers.InputLayer(input_shape=input_shape))

    # add convolutional layers
    for i in range(num_conv):
        if i==0:
            strides=[4,1]
        else:
            if padding=='same':
                strides=[4,1]
            else:
                strides=[1,1]
        model.add(
            keras.layers.Conv2D(
                channels_conv[i],
                kernel_sizes[i],
                strides,
                activation=activation,
                padding=padding,
                data_format="channels_last",
                kernel_initializer=initializer,
            )
        )
        if pooling[i] > 1:
            model.add(keras.layers.AveragePooling2D(pool_size=(1,pooling[i])))#, data_format='channels_last'))
            width_int = int((width_int - kernel_sizes[i][1] + 1)/pooling[i])
        else:
            width_int = int(width_int - kernel_sizes[i][1] + 1)
        model.add(keras.layers.Dropout(dropout_conv[i]))
    
    model.add(keras.layers.Flatten())
    
    # add fully connected layers
    for i in range(num_fc):
        model.add(
            keras.layers.Dense(
                channels_fc[i],
                activation=activation,
                kernel_initializer=initializer,
            )
        )
        model.add(keras.layers.Dropout(dropout_fc[i]))

    if task == 'mut_freq':
        model.add(
            keras.layers.Dense(1,
                activation='linear',
                kernel_initializer='glorot_normal',
            )
        )        
    if task == 'weighted_sub' or task == 'substitution':
        model.add(
            keras.layers.Dense(4,
                activation='linear',
                kernel_initializer='glorot_normal',
            )
        )        

    return model


def train_model(    # generates and trains the model with the provided parameters
    data,
    channels_conv,  #number of channels for convolutional layers (including transition layer), length N
    kernel_sizes,   #kernel sizes for convolutional layers, length N-1
    pooling,        #kernel sizes for pooling layers following each convolutional layer (pooling[i]=0 means no pooling layer there), length N-1
    dropout_conv,   #dropout rates for convolutional layers, length N-1
    channels_fc,    #number of channels for dense layers (last layer will be generated automatically based on task), length M
    dropout_fc,     #dropout rates for dense layers, length M
    task,           #can be 'mut_freq', 'substitution', 'both' or 'weighted_sub'; determines data used and architecture of the last layer
    #optimizer = tf.keras.optimizers.Adam,
    optimizer = tf.keras.optimizers.RMSprop,
    initializer = 'he_normal',
    activation = tf.nn.elu,
    save_path = False,      #where model will be saved (if not False)
    n_epochs=1,
    learning_rate=0.0001,
    mini_batch_size=32,
    loss="mse",
    log = True,
    emb = False,
): 
    [train_inputs, train_targets, test_inputs, test_targets] = data
    test_midone = give_midone(test_inputs)
    if emb == True:
        train_inputs = data_for_embedding(train_inputs)
        test_inputs = data_for_embedding(test_inputs)
        input_shape = [test_inputs.shape[1],1]
    else:
        input_shape = [test_inputs.shape[1],test_inputs.shape[2],1]
    print(input_shape)

    model = generate_model(
        input_shape,
        channels_conv,
        kernel_sizes,
        pooling, 
        dropout_conv,
        channels_fc,
        dropout_fc,
        task,
        emb = emb,
    )
    model.compile(optimizer=optimizer(learning_rate), loss=loss)
    model.fit(train_inputs, train_targets, batch_size=mini_batch_size, epochs=n_epochs)
    if save_path != False:
        model.save(save_path)
    test_loss = model.evaluate(test_inputs, test_targets, batch_size=mini_batch_size)
    
    mf_corr, sub_corr = give_corr(model.predict(test_inputs),test_targets,test_midone,task)
    return (mf_corr, sub_corr)
