#!/usr/bin/python
# @author: zhaofeng-shu33

import ace_internal
import numpy as np

def ace_cream(x, y, wt = None, delrsq = 0.01, ns = 1, cat = None):
    ''' 
    Uses the alternating conditional expectations algorithm
    to find the transformations of y and x that maximise the 
    proportion of variation in y explained by x.

    Parameters
    ----------
    x : array_like
        a matrix containing the independent variables.
        each row is an obversation of data.
    y : array_like
        a vector containing the response variable.
    wt : array_like or None, optional
        an optional vector of weights (the default is equally weighted).
    delrsq : float, optional
        termination threshold (the default is 0.01). iteration stops when rsq changes 
        less than delrsq in 3 consecutive iterations
    ns : int, optional
        number of eigensolutions (sets of transformations, the default is 1).
    cat : list
        an optional integer vector specifying which variables assume categorical values. 
        nonnegative values in cat refer to columns of the x matrix and -1 refers to the response variable.
    Returns
    -------
    tx : array_like
        the transformed x values.
    ty : array_like
        the transformed y values.

    '''
    if wt is None:
        wt = np.ones(x.shape[0])
    if(len(x.shape) == 1): # if x is rank-1 matrix
        x_internal = x.reshape([x.shape[0],1])
    else:
        x_internal = x
    x_row = x_internal.shape[0]
    x_col = x_internal.shape[1]
    iy = x_col + 1 # len(col(x)) + 1
    l = np.ones(iy,'i') # flags for each variable
    if cat is not None:
        for i in cat:
            if (i  < -1 or i >= x_col):
                raise ValueError("bad cat= specification")
            if (i == -1):
                nncol = x_col
            else:
                nncol = i
            if (l[nncol] != 5 and l[nncol] != 1):
              raise ValueError("conflicting transformation specifications")
            l[nncol] = 5
    m = np.zeros([x_row, iy],'i')
    z = np.zeros([x_row, 12], 'd')
    p = x_internal.shape[1]
    tx, ty, rsq, ierr = ace_internal.mace(x = x_internal.T, y = y, w = wt, l = l, delrsq = delrsq, m = m, z= z, ns = ns, p = p, n = x_row)
    if not(ierr == 0): # number dimension overflow
        raise ValueError("number dimension overflow")
    if(ns == 1):
        return (tx[:,:,0], ty)
    else:
        return (tx, ty)

def f_mapping(x, tx, x_candidate):
    '''
    given x and tx, find f(x) = tx and return 
    t_x_candidate = f(x_candidate)
    only used for discrete random variable X

    Parameters
    ----------
    x : array-like or list
    tx : array-like
         transformated feature of x, same dimension with x
    x_candidate : array-like or list

    Returns
    -------
    t_x_candidate : array-like

    '''
    f = {}
    feature_vector_dim = tx.shape[1]
    for i in range(len(x)):
        f[x[i]] = tx[i,:]
    t_x_candidate = np.zeros([len(x_candidate), feature_vector_dim])

    for i in range(len(x_candidate)):
        if f.get(x_candidate[i]) is None:
            t_x_candidate[i,:] = x[0,:]
        else:
            t_x_candidate[i,:] = f[x_candidate[i]]
    return t_x_candidate

if __name__ == '__main__':
    TWOPI = 2 * np.pi
    x = np.random.uniform(0,TWOPI,200)
    y = np.exp(np.sin(x)+np.random.normal(size = 200)/2)
    tx, ty = ace_cream(x,y)
    
