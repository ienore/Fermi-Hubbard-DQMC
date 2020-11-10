# -*- coding: utf-8 -*-
"""
Machine Learning based Analytical Continuation

@author: Junyin Zhang
"""
import scipy.io as sio
import numpy as np
import tensorflow as tf
matfn = u'C:/Users/ienor/Desktop/Hubbard_Advanced/MLAC/' + 'U0.mat'
batch = sio.loadmat(matfn)['mea_ML_data_p']
#mea_ML_data = NumOfVertexs*NumOfVertexs*TempSlice
#px = (px_index-NumInEdge) * delta_p - pi
NumInEdge = np.shape(batch)[1]
TempSlice = np.shape(batch)[2]

kx_index = 1
ky_index = 1
Beta = 2.0
G_k = np.zeros([1,TempSlice])
G_k[0] = batch[kx_index,ky_index,0]
for time_index in range(TempSlice-1):
    G_k[0,time_index+1] = batch[kx_index,ky_index,time_index]
    
Gauge = np.zeros([TempSlice,64])
for t_index in range(TempSlice):
    for e_index in range(64):
        tau = t_index*Beta/TempSlice
        omega = -3.2 + 0.1*e_index
        val = np.exp(-omega*(tau-Beta/2))/(2*np.cosh(Beta*omega/2))
        Gauge[t_index,e_index] = val

Gauge_tf = tf.convert_to_tensor(Gauge, dtype=tf.float32);
#################################################################
import keras
from keras import models
from keras import layers
from keras.layers import Input,Dense
from keras import backend as K


num_epochs = 30;
num_train = 1024;
num_test = 1e1;
train_data = np.zeros([round(num_train),TempSlice])
test_data = np.zeros([round(num_test),TempSlice])
for index in range(round(num_train)):
    train_data[index,:] = G_k;
for index in range(round(num_test)):
    test_data[index,:] = G_k
    
    
    

def custom_loss_wrapper(input_tensor):
    def custom_loss(y_true,y_pred):
        #G_pred = K.dot(Gauge_tf,y_pred)
        loss_value = K.sum(K.square(K.dot(Gauge_tf,K.transpose(y_pred))-K.transpose(input_tensor)))
        return loss_value
    return custom_loss

def baseline_model():
    input_tensor = Input(shape=(16,))
    hidden_1 = Dense(64, activation='sigmoid')(input_tensor)
    hidden_2 = Dense(32, activation='relu')(hidden_1)
    out = Dense(64, activation='sigmoid')(hidden_2)
    model = keras.Model(input_tensor, out)
    model.compile(optimizer='rmsprop',
                loss = custom_loss_wrapper(input_tensor))
    return model

model = baseline_model()

history = model.fit(train_data,train_data,
                    epochs = num_epochs,batch_size = 1,verbose=1)
#history = history.history['val_mean_absolute_error']


A_pred = np.transpose(model.predict(G_k))
G_pred = np.transpose(np.dot(Gauge,A_pred))
G_diff_r = (G_pred - G_k)/G_k
print(np.linalg.norm(G_diff_r))



#sample_tf = tf.convert_to_tensor(sample, dtype=tf.float32);
