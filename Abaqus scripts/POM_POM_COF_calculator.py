from scipy.interpolate import RegularGridInterpolator
import numpy as np

### This can be used to interpolate CoF from experimental results

###INPUTS:
T1 = 3
N1 = 1000

###############################
### Experimental results from:
# Walton, D., Cropper, A., Weale, D. and Klein Meuleman, P. (2002) 'The efficiency and friction of plastic cylindrical gears Part 1: Influence of materials'
torque = (3, 5, 7, 10)
rpm = (250, 500, 1000, 1500, 2000)
mu = np.zeros((4,5))
# torque  = 3
mu[0,0] = 0.38        ## 250 rpm
mu[0,1] = 0.55        ## 500 rpm
mu[0,2] = 0.67        ## 1000 rpm
mu[0,3] = 0.70        ## 1500 rpm
mu[0,4] = 0.72        ## 2000 rpm
# torque  = 5
mu[1,0] = 0.33        ## 250 rpm
mu[1,1] = 0.42        ## 500 rpm
mu[1,2] = 0.47        ## 1000 rpm
mu[1,3] = 0.48        ## 1500 rpm
mu[1,4] = 0.51        ## 2000 rpm
# torque  = 7
mu[2,0] = 0.31        ## 250 rpm
mu[2,1] = 0.35        ## 500 rpm
mu[2,2] = 0.38        ## 1000 rpm
mu[2,3] = 0.39        ## 1500 rpm
mu[2,4] = 0.40        ## 2000 rpm
# torque  = 10
mu[3,0] = 0.29        ## 250 rpm
mu[3,1] = 0.32        ## 500 rpm
mu[3,2] = 0.35        ## 1000 rpm
mu[3,3] = 0.36        ## 1500 rpm
mu[3,4] = 0.29        ## 2000 rpm

fn = RegularGridInterpolator((torque,rpm), mu)
pts = np.array([T1,N1])
print('mu = ', fn(pts))

