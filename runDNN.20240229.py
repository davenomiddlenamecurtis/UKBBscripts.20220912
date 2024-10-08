#!/usr/bin/env python3
# coding: utf-8

# run standard DNN analysis on plink files containging informative, pruned SNPs and target phenotype
# NB The test set has been used to help select SNPs so test is not independent
# To fix this, would need to either not select SNPs first or else get separate test set

# Code adapted from https://github.com/thepanlab/GattacaNet2/blob/main/Models/model_STL.py
# Thanks to adrien.f.badre-1@ou.edu for making this available

import tensorflow.compat.v1 as tf
tf.disable_v2_behavior()
import numpy as np
import sys
from pandas_plink import read_plink#1_bin
import pandas as pd
# import matplotlib.pyplot as plt # DC not using
# from sklearn.decomposition import PCA
# from sklearn.linear_model import LogisticRegression
from sklearn import metrics
from math import sqrt
import multiprocessing
from multiprocessing import Pool
import random
# from sklearn.metrics import roc_curve,roc_auc_score # DC not using

import os
# print(os.environ['CONDA_DEFAULT_ENV']) # DC does not work
# DC changed as below os.environ["CUDA_VISIBLE_DEVICES"]="1"
# https://stackoverflow.com/questions/70106418/how-can-i-run-tensorflow-without-gpu
os.environ["CUDA_VISIBLE_DEVICES"]="-1"
# DC I may not have needed this because in fact problem was low RAM

PlinkRoot="bestSNPs.HT.20240312"
(bim,fam,bed)=read_plink(PlinkRoot)
Yd=fam[['iid','trait']].astype('uint8')
# unlike the original code, I rely on all the subjects in the fam file having phenotypes so no need to get ids variable and sort genotypes with it
Yd['iid']=range(Yd.shape[0])
ids=list(set(Yd['iid']))
Yd=Yd.set_index('iid').loc[ids].values
# DC above line changes type of Yd from class 'pandas.core.frame.DataFrame' to class 'numpy.ndarray'
# and then it only has one column plus an index

Yd=Yd-1

print("Conversion")
bed=bed.astype('uint8')
bed.shape

bed=bed.T

bed.shape

x_train=bed[:bed.shape[0]*70//100,:]
x_val=bed[bed.shape[0]*70//100:bed.shape[0]*85//100,:]
x_test=bed[bed.shape[0]*85//100:,:]

y_train=Yd[:bed.shape[0]*70//100:]
assert x_train.shape[0]==y_train.shape[0]
y_val=Yd[bed.shape[0]*70//100:bed.shape[0]*85//100:]
assert x_val.shape[0]==y_val.shape[0]
y_test=Yd[bed.shape[0]*85//100::]
assert x_test.shape[0]==y_test.shape[0]


tf.compat.v1.reset_default_graph()
# X=tf.placeholder(tf.float32,shape=(None,x_train.shape[1]+1),name="X")
# DC I do not know what the +1 was for  - gave an error
# X=tf.placeholder(tf.float32,shape=(None,x_train.shape[1]),name="X")
# Y=tf.placeholder(tf.float32,shape=(None,1),name="Y")
# DC I had negative losses so I will see if this helps:
X=tf.placeholder(tf.float64,shape=(None,x_train.shape[1]),name="X")
Y=tf.placeholder(tf.float64,shape=(None,1),name="Y")


with tf.name_scope("dnn"):
    training = tf.placeholder_with_default(False, shape=(), name='training')
    initializer = tf.compat.v1.initializers.he_normal()
    hidden00_drop= tf.layers.dropout(X, 0.5, training=training)
    hidden0=tf.layers.dense(hidden00_drop, 1000, name="hidden0",activation=None, kernel_initializer=initializer)
    hidden0_norm=tf.layers.batch_normalization(hidden0, training=training, momentum=0.9)
    act_hidden0=tf.nn.leaky_relu(hidden0_norm)
    hidden0_drop = tf.layers.dropout(act_hidden0, 0.5, training=training)
    hidden1=tf.layers.dense(hidden0_drop, 200, name="hidden1",activation=None, kernel_initializer=initializer)
    hidden1_norm=tf.layers.batch_normalization(hidden1, training=training, momentum=0.9)
    act_hidden1=tf.nn.leaky_relu(hidden1_norm)
    hidden1_drop= tf.layers.dropout(act_hidden1, 0.5, training=training)
    hidden2=tf.layers.dense(hidden1_drop, 50, name="hidden2",activation=None, kernel_initializer=initializer)
    hidden2_norm=tf.layers.batch_normalization(hidden2, training=training, momentum=0.9)
    act_hidden2=tf.nn.leaky_relu(hidden2_norm)
    hidden2_drop= tf.layers.dropout(act_hidden2, 0.5, training=training)
    output=tf.layers.dense(  hidden2_drop, 1, name="output_final",activation=None)

with tf.name_scope("loss"):
    cross_entropy =tf.nn.sigmoid_cross_entropy_with_logits(labels=Y, logits=output)
    loss=tf.reduce_mean(cross_entropy)
    error=loss

with tf.name_scope("train"):
    optimizer =tf.train.AdamOptimizer(learning_rate=0.0001,beta1=0.9,beta2=0.999,epsilon=1e-08,use_locking=False,name='Adam')
    extra_update_ops = tf.get_collection(tf.GraphKeys.UPDATE_OPS)
    with tf.control_dependencies(extra_update_ops):
        training_op = optimizer.minimize(error)

with tf.name_scope("eval"):
    predicted = tf.nn.sigmoid(output)

def shuffle_batch(X, y, batch_size):
    rnd_idx = np.random.permutation(len(X))
    n_batches = len(X) // batch_size
    for batch_idx in np.array_split(rnd_idx, n_batches):
        X_batch, y_batch = X[batch_idx], y[batch_idx]
        yield X_batch, y_batch


# In[ ]:


import time
start = time.time()
batch_size=512
best_loss=7777
size_big_batch=60000
size_big_val_test=35000
saver=tf.train.Saver()
scores_test=[]
pheno_id=0

try:
    sess.close()
except:
    pass

print("Load session")
#I get a syntax error if I do not have blank line after pass
init=tf.global_variables_initializer()
loc=tf.local_variables_initializer()
# sess = tf.InteractiveSession(config=tf.ConfigProto(device_count={ "GPU":1}))
sess = tf.InteractiveSession(config=tf.ConfigProto(device_count={ "CPU": 44}))
# copied from further down because not using a GPU (yet)
init.run()
loc.run()
loss_tab=[]
epoch_tab=[]
loss_tab_val=[]

open(PlinkRoot+'_loss.txt', 'w')

for epoch in range(30):
    iteration=0
    for i in range(x_train.shape[0]//size_big_batch+1):
        x_train_big=x_train[size_big_batch*i:size_big_batch*(i+1),:].compute()
        x_train_big[np.isnan(x_train_big)]=2
        y_train_big=y_train[size_big_batch*i:size_big_batch*(i+1):]
        for x_batch,y_batch in shuffle_batch(x_train_big, y_train_big, batch_size):
            loc.run()
            sess.run(training_op,feed_dict={X:x_batch,Y:np.array([y_batch[:,pheno_id]]).T,training:True})
# got this error:  Cannot feed value of shape (513, 5523) for Tensor 'X:0', which has shape '(?, 5524)'
# before I removed "+1" from following line:
# X=tf.placeholder(tf.float32,shape=(None,x_train.shape[1]+1),name="X")
            print("%d ITERATION:%d/%d "%(epoch,iteration,len(x_train)//batch_size),end='\r')
            iteration+=1
        del x_train_big
    loc.run()
    loss_train=sess.run(loss,feed_dict={X:x_batch,Y:np.array([y_batch[:,pheno_id]]).T,training:False})
    print(epoch,"Train Loss:",loss_train)
    epoch_tab.append(epoch)
    loss_tab.append(loss_train)
    loss_values=[]
    predicted_values=[]
    for j in range(x_val.shape[0]//size_big_val_test+1):
        loc.run()
        x_val_big=x_val[size_big_val_test*j:size_big_val_test*(j+1),:].compute()
        x_val_big[np.isnan(x_val_big)]=2
        y_val_big=y_val[size_big_val_test*j:size_big_val_test*(j+1),:]
        for i in range(x_val_big.shape[0]//batch_size+1):
            x_batch_val=x_val_big[i*batch_size:(i+1)*batch_size]
            loss_val_loc,predicted_val_loc=sess.run([error,predicted],feed_dict={X:x_batch_val,Y:np.array([y_val_big[i*batch_size:(i+1)*batch_size,pheno_id]]).T,training:False})
            predicted_values.append(predicted_val_loc)
            loss_values.append(loss_val_loc)
        del x_val_big
    loss_val=np.mean(loss_values)
    pred_val_v=np.concatenate(predicted_values,axis=0)
    if loss_val<best_loss:
        save_path = saver.save(sess, PlinkRoot+"_loss.ckpt")
        best_loss=loss_val
    with open(PlinkRoot+'_loss.txt', 'a+') as file:
        file.write("Epoch"+str(epoch)+" Loss:"+str(loss_val)+"\n\n")
    print(epoch,"Loss:",loss_val)
    print("\n")
    
sess.close()
sess = tf.InteractiveSession(config=tf.ConfigProto(device_count={ "CPU": 44}))
saver.restore(sess, PlinkRoot+"_loss.ckpt")
pred_test=[]
size_big_val_test=35000
for j in range(x_test.shape[0]//size_big_val_test+1):
    loc.run()
    x_test_big=x_test[size_big_val_test*j:size_big_val_test*(j+1),:].compute()
    x_test_big[np.isnan(x_test_big)]=2
    y_test_big=y_test[size_big_val_test*j:size_big_val_test*(j+1),:]
    for i in range(x_test_big.shape[0]//batch_size+1):
        x_batch_test=x_test_big[i*batch_size:(i+1)*batch_size]
        predicted_loc=sess.run(predicted,feed_dict={X:x_batch_test,Y:np.array([y_test_big[i*batch_size:(i+1)*batch_size,0]]).T,training:False})
        pred_test.append(predicted_loc)
    del x_test_big

from sklearn import metrics
pred_test_v=np.concatenate(pred_test,axis=0)
fpr, tpr, thresholds = metrics.roc_curve(y_test[:,pheno_id],pred_test_v, pos_label=1)
precision, recall, thresholds = metrics.precision_recall_curve(y_test[:,pheno_id], pred_test_v, pos_label=1)
scores_test.append([metrics.auc(fpr, tpr),metrics.auc(recall, precision),metrics.average_precision_score(y_test[:,pheno_id],pred_test_v)])
end = time.time()
print(end - start)
   
disease='trait'

pd.DataFrame({disease:scores_test,}).to_csv(PlinkRoot+".score_test.csv")

pd.DataFrame(pred_test_v).to_csv(PlinkRoot+".prediction_test.csv")
    






