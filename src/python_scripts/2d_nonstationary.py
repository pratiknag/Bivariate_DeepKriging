#!/usr/bin/env python
# coding: utf-8

"""
Created on Tuesday August  29 23:33:04 2023

@author: Pratik
"""

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import pandas as pd
import keras
from keras.models import Sequential,Model
from keras.layers import Dense, Dropout, BatchNormalization,Input
# from keras.wrappers.scikit_learn import KerasRegressor
from keras.callbacks import EarlyStopping, ModelCheckpoint
from keras import regularizers,initializers
from keras.layers import GaussianNoise,LeakyReLU
import keras.backend as Kb
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.utils import class_weight
import numpy as np
from numpy import exp
import matplotlib.pyplot as plt
import matplotlib;matplotlib.rcParams['figure.figsize'] = (8,6)
import pylab 
import time
import sys
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import MinMaxScaler

# Load nonGaussian datasets and do classification on them 
num_sim = int(sys.argv[1])

# Functions for calculation of MSE and MAE
def mse(y_pred,y_true):
    mse = np.mean((y_pred-y_true)**2)
    return mse

def mae(y_pred,y_true):
    mae = np.mean(np.absolute(y_pred-y_true))
    return mae


def model_function(df_train, phi, sim_iteration,model_saving_location):
    
    variance_var1 = np.var(df_train["var1"])
    variance_var2 = np.var(df_train["var2"])


    def custom_mse(y_true, y_pred):

        # calculating squared difference between target and predicted values 
        loss = Kb.square(y_pred - y_true)  # (batch_size, 2)

        # multiplying the values with weights along batch dimension
        loss = loss * [(1/variance_var1), (1/variance_var2)]          # (batch_size, 2)

        # summing both loss values along batch dimension 
        loss = Kb.sum(loss, axis=1)        # (batch_size,)

        return loss

    N = len(df_train)
    df_train["var1"] = (df_train["var1"] - np.mean(df_train["var1"]))/np.sqrt(variance_var1)
    df_train["var2"] = (df_train["var2"] - np.mean(df_train["var2"]))/np.sqrt(variance_var2)
    s = np.vstack((df_train["x"],df_train["y"])).T
    y = np.array(df_train[["var1","var2"]])
    
    # DeepKriging model for continuous data
    model = Sequential()
    # model.add(Dense(100, input_dim = 2,  kernel_initializer='he_uniform', activation='relu'))
#     model.add(Dense(100, input_dim = phi.shape[1],  
#                 kernel_initializer=initializers.RandomNormal(stddev=0.01), activation='relu'))
    model.add(Dense(100, input_dim = phi.shape[1],  
                kernel_initializer='he_uniform', activation='linear'))
    model.add(LeakyReLU(alpha=0.1))
#     model.add(Dense(100, input_dim = encoder_train.shape[1],  kernel_initializer='he_uniform', activation='relu'))
    # model.add(Dropout(rate=0.5))
    # model.add(BatchNormalization())
    
#     model.add(Dense(100, kernel_regularizer=regularizers.L1L2(l1=1e-5, l2=1e-4),
#                 bias_regularizer=regularizers.l2(1e-4),
#                 activity_regularizer=regularizers.l2(1e-5),activation='relu'))
#     model.add(Dense(100, kernel_regularizer=regularizers.L1L2(l1=1e-5, l2=1e-4),
#                 bias_regularizer=regularizers.l2(1e-4),
#                 activity_regularizer=regularizers.l2(1e-5),activation='relu'))
    model.add(Dense(100, activation='linear'))
    model.add(LeakyReLU(alpha=0.1))
    model.add(Dense(100, activation='linear'))
    model.add(LeakyReLU(alpha=0.1))
    model.add(Dense(100, activation='linear'))
    model.add(LeakyReLU(alpha=0.1))
    model.add(Dense(100, activation='linear'))
    model.add(LeakyReLU(alpha=0.1))
#     model.add(Dense(100, activation='linear'))
#     model.add(LeakyReLU(alpha=0.1))
#     model.add(Dense(100, activation='relu'))
#     model.add(Dense(100, activation='relu'))
    model.add(Dense(50, activation='linear'))
    model.add(LeakyReLU(alpha=0.1))
    model.add(Dense(50, activation='linear'))
    model.add(LeakyReLU(alpha=0.1))
#     model.add(Dense(50, activation='relu'))
    # model.add(Dense(100, activation='relu'))
    #model.add(Dropout(rate=0.5))
#     model.add(Dense(50, activation='ReLU'))
#     model.add(Dense(10, activation='ReLU'))
    #model.add(BatchNormalization())
    model.add(Dense(2, activation='linear'))
    
    optimizer = keras.optimizers.Adam(learning_rate=0.001)
    model.compile(optimizer=optimizer, loss= 'mse', metrics=['mae','mse'])

    
    

    print('<<<<<<<<<<<<<<<< Fitting DNN-model for %2d-th simulation >>>>>>>>>>>>>>>>>'%(sim_iteration + 1))
    result = model.fit(phi, y, 
                       validation_split = 0.1, epochs = 500, batch_size = 256, verbose = 0)
    
    callbacks = [EarlyStopping(monitor='val_loss', patience=100),
                 ModelCheckpoint(filepath=model_saving_location+'/Biv_nonstationary_model.h5', monitor='val_loss', save_best_only=True)]
    result = model.fit(phi, y, callbacks=callbacks, 
                       validation_split = 0.1, epochs = 1000, batch_size = 256, verbose = 0)

    model = keras.models.load_model(model_saving_location+'/Biv_nonstationary_model.h5')
    return model


def main():
    print("         ")
    print("#####################################################################")
    print("############# POINT PREDICTION FOR NONSTATIONARY DATA ###############")
    print("#####################################################################")
    print("         ")
    mse_var1 = []
    mse_var2 = []
    for sim in range(num_sim):
        
        df_loc = pd.read_csv("synthetic_data_simulations_non-Stationary/2d_nonstationary_1200_"+str(sim+1)+".csv", 
                             sep = ",")
        phi = pd.read_csv("src/python_scripts/phi.csv", 
                             sep = ",")
        df_train,df_test,phi_train,phi_test = train_test_split(df_loc,pd.DataFrame.to_numpy(phi), test_size = 0.1, random_state=123)
        df_train.reset_index(drop=True, inplace=True)
        df_test.reset_index(drop=True, inplace=True)
        # Saving the training and testing datasets          
        df_train.to_csv("synthetic_data_simulations_non-Stationary/training_data/2D_nonstationary_1200_"+str(sim+1)+"-train.csv",
                                                              index = False)
        df_test.to_csv("synthetic_data_simulations_non-Stationary/testing_data/2D_nonstationary_1200_"+str(sim+1)+"-test.csv",
                        index = False)

        df_train1 = df_train.copy()
        N = len(df_train1)
        s = np.vstack((df_train1["x"],df_train1["y"])).T

        # Training the model 
        model_saving_location = 'Model_Example'
        model = model_function(df_train1,phi_train,sim,model_saving_location)

        # Basis functions for test set 

        N = len(df_test)
        s = np.vstack((df_test["x"],df_test["y"])).T
        y_test = np.array(df_test[["var1","var2"]]) 
        y_pred = model.predict(phi_test)
        y_pred[:,0] = (y_pred[:,0]*np.sqrt(np.var(df_train["var1"]))) + np.mean(df_train["var1"])
        y_pred[:,1] = (y_pred[:,1]*np.sqrt(np.var(df_train["var2"]))) + np.mean(df_train["var2"])

        mse_var1.append(mse(y_pred[:,0], y_test[:,0]))
        mse_var2.append(mse(y_pred[:,1], y_test[:,1]))
        
    df_mse = pd.DataFrame(np.vstack((mse_var1,mse_var2)).T, columns = ["mse_var1","mse_var2"])
    ## MSE saved in file 
    df_mse.to_csv("plot_results/DeepKriging_nonstationary_mse.csv", index = False)
    print(mse_var1)
    print(mse_var2)
        

if __name__ == '__main__':
    main()







