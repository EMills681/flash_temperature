from scipy.interpolate import RegularGridInterpolator
import numpy as np
import matplotlib.pyplot as plt

def interpolate_mu(T, N):
    ### Experimental results from:
    # Walton, D., Cropper, A., Weale, D. and Klein Meuleman, P. (2002) 'The efficiency and friction of plastic cylindrical gears Part 1: Influence of materials'
    torque = (3, 5, 7, 10)
    rpm = (50, 250, 500, 1000, 1500, 2000)
    mu = np.zeros((4,6))
    # torque  = 3
    mu[0,0] = 0.10        ## 50 rpm
    mu[0,1] = 0.38        ## 250 rpm
    mu[0,2] = 0.55        ## 500 rpm
    mu[0,3] = 0.67        ## 1000 rpm
    mu[0,4] = 0.70        ## 1500 rpm
    mu[0,5] = 0.72        ## 2000 rpm
    # torque  = 5
    mu[1,0] = 0.11        ## 50 rpm
    mu[1,1] = 0.33        ## 250 rpm
    mu[1,2] = 0.42        ## 500 rpm
    mu[1,3] = 0.47        ## 1000 rpm
    mu[1,4] = 0.48        ## 1500 rpm
    mu[1,5] = 0.51        ## 2000 rpm
    # torque  = 7
    mu[2,0] = 0.13        ## 50 rpm
    mu[2,1] = 0.31        ## 250 rpm
    mu[2,2] = 0.35        ## 500 rpm
    mu[2,3] = 0.38        ## 1000 rpm
    mu[2,4] = 0.39        ## 1500 rpm
    mu[2,5] = 0.40        ## 2000 rpm
    # torque  = 10
    mu[3,0] = 0.16        ## 50 rpm
    mu[3,1] = 0.29        ## 250 rpm
    mu[3,2] = 0.32        ## 500 rpm
    mu[3,3] = 0.35        ## 1000 rpm
    mu[3,4] = 0.36        ## 1500 rpm
    mu[3,5] = 0.29        ## 2000 rpm

    ############################################
    fn = RegularGridInterpolator((torque,rpm), mu)
    pts = np.array([T,N])

    return fn(pts)

print(interpolate_mu(3, 1000))
T1 = 3
T2 = 5
T3 = 7
T4 = 10
T5 = 4
T6 = 6
T7 = 8
T8 = 9

steps = 11

N1 = np.linspace(50,2000, steps)
mu1 = np.zeros(steps)
mu2 = np.zeros(steps)
mu3 = np.zeros(steps)
mu4 = np.zeros(steps)

mu5 = np.zeros(steps)
mu6 = np.zeros(steps)
mu7 = np.zeros(steps)
mu8 = np.zeros(steps)

for i in range(steps):
    mu1[i] = interpolate_mu(T1, N1[i])
    mu2[i] = interpolate_mu(T2, N1[i])
    mu3[i] = interpolate_mu(T3, N1[i])
    mu4[i] = interpolate_mu(T4, N1[i])

    mu5[i] = interpolate_mu(T5, N1[i])
    mu6[i] = interpolate_mu(T6, N1[i])
    mu7[i] = interpolate_mu(T7, N1[i])
    mu8[i] = interpolate_mu(T8, N1[i])

plt.plot(N1,mu1, marker = 'D', label="3 N m")
plt.plot(N1,mu5, '--', marker = 'o', label="4 N m")
plt.plot(N1,mu2, marker = 's', label="5 N m")
plt.plot(N1,mu6, '--', marker = 'p', label="6 N m")
plt.plot(N1,mu3, marker = '^', label="7 N m")
plt.plot(N1,mu7, '--', marker = '*', label="8 N m")
plt.plot(N1,mu8, '--', marker = 'v', label="9 N m")
plt.plot(N1,mu4, marker = 'x', label="10 N m")

plt.xlim(0,2500)
plt.ylim(bottom = 0)
plt.legend(loc="lower right")
plt.xlabel('Rotational Speed (rpm)')
plt.ylabel('Coefficient of Friction')
plt.grid()
plt.rc('font', size = 20)
plt.show()
