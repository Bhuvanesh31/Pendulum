import numpy as np
import matplotlib.pyplot as plt


theta_c = 1
theta_dot_c = 0.1
dt = 0.01
steps = 1000
g = 9.81
l_1 = 1
time = np.arange(0,dt*steps,dt)
alphas = np.arange(0,1,1/steps)
rgba_colors = np.zeros((100,4))
rgba_colors[:,0] = 0.2
rgba_colors[:,1] = 0.7
rgba_colors[:,2] = 0.5
rgba_colors[:,3] = np.linspace(0.,1.,100)
theta = np.zeros(steps)
theta_dot = np.zeros(steps)

plt.plot(rgba_colors[-1000:0:-10,3])
plt.show()


plt.axis([0, steps*dt, -np.pi, np.pi])
plt.ion()
plt.show()

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
    return result

for i in range(steps):
    theta[i], theta_dot[i] = exp_Eul((DGL_single_theta, DGL_single_theta_dot), (theta_c, theta_dot_c))
    theta_c = theta[i]
    theta_dot_c = theta_dot[i]

    if not i%10:
        a = (i>100)*(i-100)
        b = len(time[a:i])-1
        print(len(time[a:i]))
        print(a, i)
        if i > 0:
            pts.remove()
        pts = plt.scatter(time[a:i:10], theta[a:i:10], color=rgba_colors[-b::10,:])


        plt.pause(0.0001)
