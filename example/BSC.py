#!/usr/bin/python
#author: zhaofeng-shu33

import numpy as np
from ace_cream import ace_cream

def pearson_correlation(X,Y):
    return (np.mean(X*Y, axis=0) -np.mean(X, axis = 0)* np.mean(Y, axis = 0)) / ( np.std(X, axis = 0) * np.std(Y, axis = 0))

if __name__ == '__main__':
    N_SIZE = 1000
    P_CROSSOVER = 0.1
    x = np.random.choice([0,1],size=N_SIZE)
    y = np.random.uniform(size=N_SIZE)
    y = y > (1- P_CROSSOVER)
    y = y*(1-x) + (1-y)*x
    # x = a_1 * x + a_2
    # y = b_1 * y + b_2 
    # a_1, a_2, b_1, and b_2 are constants has no influence on the result
    # by the linearity of pearson correlation
    # We conclude that the alphabet is irrelevant if the size is 2
    # theoretical result = 2 * max(p, 1-p) - 1
    print('rho(x,y)',pearson_correlation(x,y))
    # use fortran ace by 1985 article author
    tx, ty = ace_cream(x, y, cat = [-1,0])
    print('mean(tx) = %f, std(tx) = %f'%(np.mean(tx), np.std(tx)))
    print('mean(ty) = %f, std(ty) = %f'%(np.mean(ty), np.std(ty)))
    print('rho(tx,ty)',pearson_correlation(tx,ty))
    # as can be seen, the transformation is simply a normalization trick
    # tx = \frac{x-E[X]}{\sqrt{\Var[X]}}

