# -*- coding: utf-8 -*-
"""
Reinforcement Learning Enhanced Analytic Continuation
RLEAC

@author: Junyin Zhang (Discussed with Bo Zhang)
"""

import numpy as np
from matplotlib import pyplot as plt
TempSlice = 16
Beta = 2.0
D_Tau = Beta / TempSlice
N_epoch = int(1e1)
N_warm = int(1e5)
step = 1e-2
E_min = -20
E_max = 20
dE = 0.2
E_range = np.arange(E_min,E_max+dE,dE)
G_k = np.zeros([TempSlice,1]);
Gauge = np.zeros([TempSlice,len(E_range)])
Theta = 1e6;
diff = 1.0
#######################################################
A = np.zeros([len(E_range),1])
for index in range(len(E_range)):
    omega = E_range[index]
    A[index] = np.exp(-0.2*omega**2)
    
A = A/(np.sum(A)*dE)
D_Tau = Beta/TempSlice
A_ori = A*1.0
for time_index in range(TempSlice):
    for E_index in range(len(E_range)):
        tau = (time_index-1)*D_Tau
        omega = E_range[E_index]
        Gauge[time_index,E_index] = dE * np.exp(-omega*(tau - Beta/2))/(2*np.cosh(Beta*omega/2))

G_k = np.dot(Gauge,A)
for E_index in range(len(E_range)):
    A[E_index] = np.random.rand()

A = A/(np.sum(A)*dE)
############################## Calculation Process ############################3
A_collect = np.zeros([len(E_range),N_epoch])


##################################################################################
import keras
import tensorflow as tf
from keras import models
from keras import layers
from keras.layers import Input,Dense
from keras import backend as K

num_epochs = 1;
num_train = 1024;
num_test = 1e1;
Gauge_tf = tf.convert_to_tensor(Gauge, dtype=tf.float32);


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
    out = Dense(len(E_range), activation='sigmoid')(hidden_2)
    model = keras.Model(input_tensor, out)
    model.compile(optimizer='rmsprop',
                loss = custom_loss_wrapper(input_tensor))
    return model


model = baseline_model()

train_data = np.zeros([round(num_train),TempSlice])
test_data = np.zeros([round(num_test),TempSlice])
for epoch_index in range(N_epoch):
    model.fit(train_data,train_data,
              epochs = num_epochs,batch_size = 1,verbose=1)
    A_new = np.transpose(model.predict(G_k.reshape([1,16])))
    value_old = np.sum((np.dot(Gauge,A) - G_k)**2)
    value_new = np.sum((np.dot(Gauge,A_new) - G_k)**2)
    diff = np.linalg.norm(np.dot(Gauge,A) - G_k)
    step = np.sqrt(diff)
    if np.random.rand()<np.exp(Theta*(value_old - value_new)) and np.min(A_new) >=0:
        A = A_new*1.0
    A_collect[:,epoch_index] = A.reshape(len(E_range))





#plt.plot(E_range,A_final)
#plt.plot(E_range,A_ori)