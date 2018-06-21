#!/usr/bin/python
#author: zhaofeng-shu33

import numpy as np
from ace_cream import ace_cream

def pearson_correlation(X,Y):
    return (np.mean(X*Y, axis=0) -np.mean(X, axis = 0)* np.mean(Y, axis = 0)) / ( np.std(X, axis = 0) * np.std(Y, axis = 0))

if __name__ == '__main__':
    N_SIZE = 1000
    ERROR_PROBABILITY = 0.1
    x = np.random.choice([0,1],size=N_SIZE)
    y = np.random.uniform(size=N_SIZE)
    for i in range(len(x)):
        if(y[i] < ERROR_PROBABILITY):
            y[i] = 2
        else:
            y[i] = x[i]
    dic_Y = {0:6, 1:8, 2:3}
    dic_X = {0:7, 1:9}
    for i in range(len(y)):
        y[i] = dic_Y[y[i]]
        x[i] = dic_X[x[i]]
    print('rho(x,y)',pearson_correlation(x,y))
    # use fortran ace by 1985 article author
    tx, ty = ace_cream(x, y, cat = [-1,0])
    print('mapped X symbol list: ')
    print(np.unique(tx))
    print('mapped Y symbol list: ')
    print(np.unique(ty))

    print('mean(tx) = %f, std(tx) = %f'%(np.mean(tx), np.std(tx)))
    print('mean(ty) = %f, std(ty) = %f'%(np.mean(ty), np.std(ty)))
    print('rho(tx,ty)',pearson_correlation(tx,ty))

