import numpy as np
import matplotlib.pyplot as plt


theta_c = 1
theta_dot_c = 0.1
dt = 0.001
steps = 10000
g = 9.81
l_1 = 1

theta = np.zeros(steps)
theta_dot = np.zeros(steps)

def DGL_single_theta(params):
    theta_dot = params[1]

    return  theta_dot

def DGL_single_theta_dot(params):
    theta = params[0]
    return -g/l_1 * np.sin(theta)

def exp_Eul(DGLs, vals):
    result = []
    for i, DGL in enumerate(DGLs):
        result.append(vals[i] + dt*DGL(vals))
    print(result)
    return result

for i in range(steps):
    theta[i], theta_dot[i] = exp_Eul((DGL_single_theta, DGL_single_theta_dot), (theta_c, theta_dot_c))

    theta_c = theta[i]
    theta_dot_c = theta_dot[i]

plt.plot(theta)
plt.show()
