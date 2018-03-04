#----------- Imports ------------

import numpy as np
import matplotlib.pyplot as plt
import time as tm

#----------- Constants and stuff ------------

dt = 0.0001 # size of time steps
t_max = 10  # maximum time
steps = int(t_max/dt)   # number of steps
g = 9.81    # gravitational acceleration
l_1 = 1     # length of pendulum
interval = .1/dt    # plot update interval
plot_length = int(1/dt) # number of time steps between first and last point in plot
plot_points = 20    # number of points in plot
point_interval = plot_length//plot_points   # step interval between points in plot

 #----------- Initial Conditions ------------

theta_c = np.pi/2+1 # IC for theta
theta_dot_c = 0.0   # IC for theta_dot

#----------- Arrays ------------

time = np.arange(0,dt*steps,dt) # array of time points
rgba_colors = np.zeros((plot_length,4)) # array for color used in plot
rgba_colors[:,0] = 0.2
rgba_colors[:,1] = 0.6
rgba_colors[:,2] = 0.3
rgba_colors[:,3] = np.linspace(0.,1.,plot_length)
theta = np.zeros(steps) # theta at every time step
theta_dot = np.zeros(steps) # theta_dot at every time step
x = np.zeros(steps) # x-coordinate for every theta
y = np.zeros(steps) # y-coordinate for every theta

#----------- Functions ------------

def DE_single_theta(params):
    #Differential equation for theta in a single pendulum
    theta_dot = params[1]
    return  theta_dot

def DE_single_theta_dot(params):
    #Differential equation for theta_dot in a single pendulum
    theta = params[0]
    return -g/l_1 * np.sin(theta)

def exp_Eul(DAEs, vals):
    ''' Explicit-Euler solver for system of differential equations.
        Arguments:
            DAEs    tuple of Differential Equations
            vals    tuple of current values
    '''
    result = []
    for i, DE in enumerate(DAEs):
        result.append(vals[i] + dt*DE(vals))
    return result

def RK4(DAEs, vals):
    ''' Runge-Kutta-4 solver for system of differential equations.
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

#plt.axis([0, steps*dt, -np.pi, np.pi]) # for theta / t plot
plt.axis([-1.1*l_1, 1.1*l_1, -1.1*l_1, 1.1*l_1])    # plot limits
plt.ion()   # interactive plot
plt.show()

timer_start = tm.time() # timer-start
for i in range(steps):

    theta[i], theta_dot[i] = RK4((DE_single_theta, DE_single_theta_dot), (theta_c, theta_dot_c))    # calculate next values
    theta_c = theta[i]  # overwrite old values
    theta_dot_c = theta_dot[i]
    x[i], y[i] = -l_1*np.cos(theta[i]), l_1*np.sin(theta[i])    # calculate x and y from generalized coordinate theta

    if not i%interval:  # plot every [interval] steps
        start = (i>plot_length)*(i-plot_length) # earliest point shown in plot
        if i > 0:   # remove old points and line
            pts.remove()
            bar.pop(0).remove()
        #pts = plt.scatter(time[start:i:point_interval], theta[start:i:point_interval], color=rgba_colors[-(i-start)::point_interval,:])    # for theta / t plot
        pts = plt.scatter(y[start:i:point_interval], x[start:i:point_interval], color=rgba_colors[-(i-start)::point_interval,:])    # plot points
        bar = plt.plot((0,y[i-point_interval+1]),(0,x[i-point_interval+1]), color=rgba_colors[plot_length-1,:]) # plot pendulum bar
        plt.scatter([0],[0], color=rgba_colors[-1,:])   # draw pendulum mount
        plt.pause(0.0000001)

print("duration: ", tm.time()-timer_start)  # print duration of program run
