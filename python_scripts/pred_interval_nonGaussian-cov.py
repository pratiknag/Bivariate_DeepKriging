#!/usr/bin/env python
# coding: utf-8

"""
Created on Tuesday August  29 23:33:04 2023

@author: Pratik
"""

import pandas as pd
import keras
from keras.models import Sequential,Model
from keras.layers import Dense, Dropout, BatchNormalization,Input
from keras.wrappers.scikit_learn import KerasRegressor
from keras.callbacks import EarlyStopping, ModelCheckpoint
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import Pipeline
import numpy as np
from numpy import exp
# Library for Gaussian process
# import GPy
##Library for visualization
import matplotlib.pyplot as plt
# %matplotlib inline
# %config InlineBackend.figure_format = 'svg'
# import matplotlib;matplotlib.rcParams['figure.figsize'] = (8,6)
from sklearn.model_selection import train_test_split
import pylab 
import time
from tqdm import tqdm
from python_libs.pred_interval_functions import mse,mae,fit_model,fit_ensemble,y_list_uni,variance,predict_with_pi
from python_libs.pred_interval_functions import calc_distance,get_nearest_data


# number of simulations
num_sim = 30

def main():
    print("         ")
    print("#####################################################################")
    print("####### PREDICTION INTERVALS FOR NONGAUSSIAN WITH COVARIATES ########")
    print("#####################################################################")
    print("         ")
    mpiw_var1 = []
    mpiw_var2 = []
    picp_var1 = []
    picp_var2 = []
    for sim in range(num_sim):
        df_train = pd.read_csv("synthetic_data_simulations_nonGaussian_cov/training_data/2D_nonGaussian_1200_"+str(sim+1)+"-train.csv", 
                             sep = ",")
        df_test = pd.read_csv("synthetic_data_simulations_nonGaussian_cov/testing_data/2D_nonGaussian_1200_"+str(sim+1)+"-test.csv", 
                             sep = ",")
        
        N_train = len(df_train)
        s_train = np.vstack((df_train["x"],df_train["y"])).T
        variance_var1 = np.var(df_train["var1"])
        variance_var2 = np.var(df_train["var2"])
        mean1 = np.mean(df_train["var1"])
        mean2 = np.mean(df_train["var2"])
        df_train["var1"] = (df_train["var1"] - mean1)/np.sqrt(variance_var1)
        df_train["var2"] = (df_train["var2"] - mean2)/np.sqrt(variance_var2)
        y_train = np.array(np.array(df_train[["var1","var2"]]))
        covariates_train = np.array(df_train[["cov1","cov2",
                                             "cov3","cov4",
                                             "cov5"]])
        scaler = MinMaxScaler()
        scaler.fit(covariates_train)
        covariates_train = scaler.transform(covariates_train)
        
        N_test = len(df_test)
        s_test = np.vstack((df_test["x"],df_test["y"])).T
        y_test = np.array(np.array(df_test[["var1","var2"]]))
        
        num_basis = [3**2,7**2,11**2]
        knots_1d = [np.linspace(0,1,int(np.sqrt(i))) for i in num_basis]
        ##Wendland kernel
        K = 0
        phi = np.zeros((N_train, sum(num_basis)))

        for res in range(len(num_basis)):
            theta = 1/np.sqrt(num_basis[res])*2.5
            knots_s1, knots_s2 = np.meshgrid(knots_1d[res],knots_1d[res])
            knots = np.column_stack((knots_s1.flatten(),knots_s2.flatten()))
            for i in range(num_basis[res]):
                d = np.linalg.norm(s_train-knots[i,:],axis=1)/theta
                for j in range(len(d)):
                    if d[j] >= 0 and d[j] <= 1:
                        phi[j,i + K] = (1-d[j])**6 * (35 * d[j]**2 + 18 * d[j] + 3)/3
                    else:
                        phi[j,i + K] = 0
            K = K + num_basis[res]
        phi = np.hstack((covariates_train,phi))
        
        K = 0
        phi_test = np.zeros((N_test, sum(num_basis)))

        for res in range(len(num_basis)):
            theta = 1/np.sqrt(num_basis[res])*2.5
            knots_s1, knots_s2 = np.meshgrid(knots_1d[res],knots_1d[res])
            knots = np.column_stack((knots_s1.flatten(),knots_s2.flatten()))
            for i in range(num_basis[res]):
                d = np.linalg.norm(s_test-knots[i,:],axis=1)/theta
                for j in range(len(d)):
                    if d[j] >= 0 and d[j] <= 1:
                        phi_test[j,i + K] = (1-d[j])**6 * (35 * d[j]**2 + 18 * d[j] + 3)/3
                    else:
                        phi_test[j,i + K] = 0
            K = K + num_basis[res]
        covariates_test = np.array(df_test[["cov1","cov2",
                                             "cov3","cov4",
                                             "cov5"]])
        covariates_test = scaler.transform(covariates_test)
        phi_test = np.hstack((covariates_test,phi_test))
        
        s_train_ensemble, s_train_mse, X_train_ensemble, X_train_mse, y_train_ensemble, y_train_mse= train_test_split(s_train,         phi, y_train, test_size=0.1)
        
        # base model layers 

        input_dim = Input(shape = (phi.shape[1], ))
        layer1 = Dense(100, kernel_initializer='he_uniform', activation = 'relu')(input_dim)
        layer2 = Dense(100, activation = 'relu')(layer1)
        layer3 = Dense(100, activation = 'relu')(layer2)
        layer4 = Dense(50, activation = 'relu')(layer3)
        layer5 = Dense(50, activation = 'relu')(layer4)
        final_layer = Dense(2, activation = 'linear')(layer5)
        
        s_train_ensemble1, _, X_train_ensemble1, __, y_train_ensemble1, ___= train_test_split(s_train_ensemble, X_train_ensemble, y_train_ensemble, test_size=0.2)
        
        base_model = Model(inputs = input_dim, outputs = final_layer)
        optimizer = keras.optimizers.Adam(lr=0.01)
        # Compile the Model
        base_model.compile(optimizer = optimizer, loss = 'mae')
#         base_model.summary()
        print('<<<<<<<<<<<<<<<< Fitting Base-model for %2d-th simulation >>>>>>>>>>>>>>>>>'%(sim + 1))
        base_model.fit(X_train_ensemble1, y_train_ensemble1,validation_split = 0.1, epochs = 500,
                        batch_size = 64,verbose = 0)
        base_model1 = Model(inputs = input_dim, outputs = layer3)
        
        # fit ensemble
        # start_time = time.time()
        n_members = 20
        print('********************* Fitting ensamble-models *********************')
        ensemble = fit_ensemble(n_members, X_train_ensemble, y_train_ensemble, base_model1)
        print(' --- Finished ---')
        
        # train data mean and variance vectors
        mean_vec1, var_vec1, mean_vec2, var_vec2 = predict_with_pi(ensemble, X_train_mse)
        
        # random error calculation on training data
        r1 = list()
        r2 = list()
        r1 = (y_train_mse[:,0] - mean_vec1)**2 - var_vec1
        for i in range(len(r1)):
            if r1[i] < 0: r1[i] = 0.0
        r2 = (y_train_mse[:,1] - mean_vec2)**2 - var_vec2
        for i in range(len(r2)):
            if r2[i] < 0: r2[i] = 0.0
                
        # random error calculation on test data with neighbourhood approach
        print('********************* Calculating variance of nuggets *********************')
        r1_pred = get_nearest_data(s_train_mse,s_test,r1,40)
        r2_pred = get_nearest_data(s_train_mse,s_test,r2,40)
        print(' --- Finished ---')
        
        # mean and variance vector for prediction data
        mean_vec1, var_vec1, mean_vec2, var_vec2 = predict_with_pi(ensemble, phi_test)

        var_vec1 = np.asarray(var_vec1)
        var_vec2 = np.asarray(var_vec2)
        # e = np.asarray(e)
        lower_bound1 = (np.asarray(mean_vec1)*np.sqrt(variance_var1) + mean1) - 1.96*np.asarray(np.sqrt(var_vec1 + r1_pred))*np.sqrt(variance_var1)
        upper_bound1 = (np.asarray(mean_vec1)*np.sqrt(variance_var1) + mean1) + 1.96*np.asarray(np.sqrt(var_vec1 + r1_pred))*np.sqrt(variance_var1)
        lower_bound2 = (np.asarray(mean_vec2)*np.sqrt(variance_var2) + mean2) - 1.96*np.asarray(np.sqrt(var_vec2 + r2_pred))*np.sqrt(variance_var2)
        upper_bound2 = (np.asarray(mean_vec2)*np.sqrt(variance_var2) + mean2) + 1.96*np.asarray(np.sqrt(var_vec2 + r2_pred))*np.sqrt(variance_var2)
        
        count_var0 = 0
        count_var1 = 0
        for i in range(len(y_test)):
            if ((y_test[i,0] > lower_bound1[i]) and (y_test[i,0] < upper_bound1[i])): count_var0 +=1
            if y_test[i,1] > lower_bound2[i] and y_test[i,1] < upper_bound2[i]: count_var1 +=1
        picp1 = count_var0/len(y_test)
        picp2 = count_var1/len(y_test)
        print('PICP variable 1 %4f, variable 2 %4f'%(picp1,picp2))
        width1 = np.mean(upper_bound1 - lower_bound1)
        width2 = np.mean(upper_bound2 - lower_bound2)
        print('MPIW variable 1 %4f, variable 2 %4f'%(width1,width2))
        mpiw_var1.append(width1)
        mpiw_var2.append(width2)
        picp_var1.append(picp1)
        picp_var2.append(picp2)
        
    df_mpiw = pd.DataFrame(np.vstack((mpiw_var1,mpiw_var2,picp_var1,picp_var2)).T, columns = ["mpiw_var1","mpiw_var2","picp_var1","picp_var2"])
    df_mpiw.to_csv("plot_results/DeepKriging_nongaussian_interval.csv", index = False)
    print('Average PICP variable 1 %4f, variable 2 %4f'%(np.mean(picp_var1),np.mean(picp_var2)))
    print('Average MPIW variable 1 %4f, variable 2 %4f'%(np.mean(mpiw_var1),np.mean(mpiw_var2)))
    
if __name__ == '__main__':
    main()
    








