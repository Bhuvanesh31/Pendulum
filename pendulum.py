#----------- Imports ------------

import numpy as np
import matplotlib.pyplot as plt
import time as tm

#----------- Constants and Initial Values ------------

theta_c = np.pi/2+1
theta_dot_c = 0.0

dt = 0.0001
t_max = 10
steps = int(t_max/dt)
g = 9.81
l_1 = 1
interval = .1/dt
plot_length = int(1/dt)
plot_points = 20
point_interval = plot_length//plot_points

#----------- Arrays ------------

time = np.arange(0,dt*steps,dt)
alphas = np.arange(0,1,1/steps)
rgba_colors = np.zeros((plot_length,4))
rgba_colors[:,0] = 0.2
rgba_colors[:,1] = 0.6
rgba_colors[:,2] = 0.3
rgba_colors[:,3] = np.linspace(0.,1.,plot_length)
theta = np.zeros(steps)
theta_dot = np.zeros(steps)
x = np.zeros(steps)
y = np.zeros(steps)

#----------- Functions ------------

def DE_single_theta(params):
    theta_dot = params[1]
    return  theta_dot

def DE_single_theta_dot(params):
    theta = params[0]
    return -g/l_1 * np.sin(theta)

def exp_Eul(DAEs, vals):
    result = []
    for i, DE in enumerate(DAEs):
        result.append(vals[i] + dt*DE(vals))
    return result

def RK4(DAEs, vals):
    ''' Runge-Kutta-4 solver for Ssstem of differential equations.
        Arguments:
            DAEs    tuple of Differential Equations
            vals    tuple of current values
    '''
    k = np.zeros((4,len(DAEs)))

    k[0] = [DE(vals) for DE in DAEs]
    k[1] = [DE(vals+dt*k[0]/2) for DE in DAEs]
    k[2] = [DE(vals+dt*k[1]/2) for DE in DAEs]
    k[3] = [DE(vals+dt*k[2]) for DE in DAEs]

    results = vals + dt/6*(k[0]+2*k[1]+2*k[2]+k[3])
    return tuple(results)
#----------- Execution ------------

#plt.axis([0, steps*dt, -np.pi, np.pi])
plt.axis([-1.1*l_1, 1.1*l_1, -1.1*l_1, 1.1*l_1])
plt.ion()
plt.show()

asd = tm.time()
for i in range(steps):

    theta[i], theta_dot[i] = RK4((DE_single_theta, DE_single_theta_dot), (theta_c, theta_dot_c))
    theta_c = theta[i]
    theta_dot_c = theta_dot[i]
    x[i], y[i] = -l_1*np.cos(theta[i]), l_1*np.sin(theta[i])
'''
    if not i%interval:
        start = (i>plot_length)*(i-plot_length)

        if i > 0:
            pts.remove()
            bar.pop(0).remove()
        #pts = plt.scatter(time[start:i:point_interval], theta[start:i:point_interval], color=rgba_colors[-(i-start)::point_interval,:])
        pts = plt.scatter(y[start:i:point_interval], x[start:i:point_interval], color=rgba_colors[-(i-start)::point_interval,:])
        bar = plt.plot((0,y[i-point_interval+1]),(0,x[i-point_interval+1]), color=rgba_colors[plot_length-1,:])
        plt.scatter([0],[0], color=rgba_colors[-1,:])
        plt.pause(0.0001)
'''
print("duration: ", tm.time()-asd)
