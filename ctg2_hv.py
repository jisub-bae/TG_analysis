import os
import numpy as np
import pandas as pd
from pandas import DataFrame as df
import matplotlib
from matplotlib import pyplot as plt
from scipy.io import wavfile
from scipy import signal
import scipy.io
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis as LDA
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis as QDA
from sklearn import svm
from sklearn.linear_model import LogisticRegression as LG
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import KFold

from sklearn.metrics import roc_curve,roc_auc_score
from sklearn.metrics import auc

import itertools
import pickle
import joblib
import random

Ts=500 #sampling rate (hz)
window=20 #window duration (ms), need to condider sampling rate
shift_w=2 #moving window (ms), need to condider sampling rate

window=int((window*Ts)/1000)
shift_w=int((shift_w*Ts)/1000)
if window<1:
    window=1
if shift_w<1:
    shift_w=1


def mdd(Xit2,Yit2,proto_clf,t):
     
  #  tt=len(DB[19][4][0,:,0])
 #   pp=len(DB)
#    kk=len(DB[0])

    tt=750
    
    a3_out=[]
    for tf in t:
        
        clf=pickle.loads(proto_clf[tf])

        a2_out=[]
        for t2 in range(0,tt-window,shift_w):

            X=Xit2[t2]
            y=Yit2[t2]

            pout=clf.predict(X)
            a_out=roc_auc_score(y, pout)

            a2_out.append(a_out)

        a3_out.append(a2_out)

    return a3_out
