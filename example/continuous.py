#!/usr/bin/python
# author: zhaofeng-shu33
# file-description: show ace_cream power in continuous case
import numpy as np
from ace_cream import ace_cream
HAS_MATPLOTLIB = True
try:
    import matplotlib.pyplot as plt
except Exception as e:
    HAS_MATPLOTLIB = False
    pass
if __name__ == '__main__':
    x = np.random.uniform(0,2*np.pi,400)
    y = np.exp(np.sin(x)+np.random.normal(size = 400)/2)
    print(np.corrcoef(x,y)[0,1])
    tx, ty = ace_cream(x, y)
    print(np.corrcoef(tx[:,0],ty[:,0])[0,1])
    if(HAS_MATPLOTLIB):
        plt.subplot(2,1,1)
        plt.title('np.exp(np.sin(x))')
        plt.scatter(x,tx[:,0])
        plt.xlabel('x')
        plt.ylabel('tx')
        plt.subplot(2,1,2)
        plt.scatter(y,ty[:,0])
        plt.xlabel('y')
        plt.ylabel('ty')
        plt.savefig('continuous.svg')
