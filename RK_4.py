#Importing necessary modules
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jn, yn

#Function to return z and dz/dx to use in rk4
def f(T0, m, x, y, n, L,w):
    gamma = 1.4
    R = 287
    T = T0 + m * x
    p, z = y  
    dz_dx = -(1 / T) * m * z - p * (w**2/(gamma*R*T))
    return [z, dz_dx]  

#RK4 method
def rk_4(f, x, y0, h, m, T0, L, n, w):
    p, z = [y0[0]], [y0[1]]
   
    for i in range(int(L / h)):
        p_i, z_i  = p[i], z[i] 
        k1 =  h*np.array(f(T0, m, x[i], [p_i, z_i], n, L,w))
        k2 =  h*np.array(f(T0, m, x[i] + h / 2, [p_i + k1[0] / 2, z_i + k1[1] / 2], n, L,w))
        k3 =  h*np.array(f(T0, m, x[i] + h / 2, [p_i + k2[0] / 2, z_i + k2[1] / 2], n, L,w))
        k4 =  h*np.array(f(T0, m, x[i] + h, [p_i + k3[0], z_i + k3[1]], n, L,w))
        p_i_1 = p_i + (1 / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        z_i_1 = z_i + (1 / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        p.append(p_i_1)
        z.append(z_i_1)

    return p, z

#Function to use RK4 for each omega
def check(f,x,y0,h,m,T0,L,n,w):
    p , z = [] , []

    for i in range(len(x)):
        p_i , z_i = rk_4(f, x, y0, h, m, T0, L, n, w[i])
        p.append(p_i)
        z.append(z_i)
    return p, z

P0 = 101325 #Using Patm as the mean pressure

#######  PART-A of code for finding frequencies and plotting p and u for all T0 and m by varying n

p_all = []
u_all = []

T0 = [300,500,700,900,1100]
m  = [0,-50,-100,-150,-200]

#Looping through T0 and m to get frequencies and plots for all conditions
for j in range(5):
    n = 1   # Keep changing n. n = 1,3,5.....
    L = 4
    h = 0.1
    gamma = 1.4
    R = 287
    x = np.arange(0,L+h,h)
    y0 = [2000,0]
    T = []
    w = []

#Computing Temperature and omegas at each location
    for i in range(len(x)):
        T.append(T0[j] + m[j]*x[i])
        w.append(math.pi * n * math.sqrt(gamma*R*T[i])/(2 * L))

    pressure , z = check(f,x,y0,h,m[j],T0[j],L,n,w)
    min_P = pressure[0]
    iter = []

#Simple sorting to find out which omega gives minimum pressure at x = L
    for i in range(len(x)):
        if abs(pressure[i][-1]) <= abs(min_P[-1]):
            min_P = pressure[i]
            iter.append(i)
  
    z_min = z[iter[-1]]
    x_min = x[iter[-1]]

    T_min = T0[j] + m[j]*x_min
    f_min = n*math.sqrt(gamma*R*T_min)/(4*L)
    print("frequency: ",f_min)

    u = []
    d = []

# Finding velocity from momentum equation
    for i in range(len(x)):
        d_temp = (P0/(R*T[i]))
        u.append((-1*z_min[i])/(1j * 2 * math.pi * f_min * d_temp))   
        d.append(d_temp)

#Finding absolute values of pressure and velocity
    p_abs =  [abs(x) for x in min_P]
    p_all.append(p_abs)
    u_abs =  [abs(x) for x in u]
    u_all.append(u_abs)

for i in range(5):
    plt.plot(x,p_all[i],label=f"T0 = {T0[i]}, m = {m[i]} n = 1")
    plt.xlabel("Distance (m)")
    plt.ylabel("Acoustic Pressure Amplitude (Pa)")
    plt.title("Acoustic Pressure Amplitude vs x")
    plt.legend(loc = 'best', fontsize="small", prop={'size': 6})
    plt.grid()
plt.show()

for i in range(5):
    plt.plot(x,u_all[i],label=f"T0 = {T0[i]}, m = {m[i]} n = 1")
    plt.xlabel("Distance (m)")
    plt.ylabel("Acoustic Velocity Amplitude (m/s)")
    plt.title("Acoustic Velocity Amplitude vs x")
    plt.legend(loc = 'best', fontsize="small", prop={'size': 6})
    plt.grid()
plt.show()

###### Part-B of code for reproducing the graphs for two separate initial conditions ########

T02 = 1100
m2 = -200
T01 = 500
m1 = -50 
n_1 = 5
n_2 = 5 
L = 4
h = 0.01
gamma = 1.4
R = 287
x = np.arange(0,L+h,h)
y0 = [2000,0]
T1 = []
T2 = []
w_1 = []
w_2 = []

# Finding T and omega at each location
for i in range(len(x)):
    T1.append(T01 + m1*x[i])
    T2.append(T02 + m2*x[i])
    w_1.append(math.pi * n_1 * math.sqrt(gamma*R*T1[i])/(2 * L))
    w_2.append(math.pi * n_2 * math.sqrt(gamma*R*T2[i])/(2 * L))


pressure_1 , z_1 = check(f,x,y0,h,m1,T01,L,n_1,w_1)
min_P_1 = pressure_1[0]
iter_1 = []

#Simple sorting to find out which omega gives minimum pressure at x = L
for i in range(len(x)):
    if abs(pressure_1[i][-1]) < abs(min_P_1[-1]):
        min_P_1 = pressure_1[i]
        iter_1.append(i)

z_min_1 = z_1[iter_1[-1]]
x_min_1 = x[iter_1[-1]]

T_min_1 = T01 + m1*x_min_1
f_min_1 = n_1*math.sqrt(gamma*R*T_min_1)/(4*L)

pressure_2 , z_2 = check(f,x,y0,h,m2,T02,L,n_2,w_2)
min_P_2 = pressure_2[0]
iter_2 = []

#Simple sorting to find out which omega gives minimum pressure at x = L
for i in range(len(x)):
    if abs(pressure_2[i][-1]) <= abs(min_P_2[-1]):
        min_P_2 = pressure_2[i]
        iter_2.append(i)

z_min_2 = z_2[iter_2[-1]]
x_min_2 = x[iter_2[-1]]

T_min_2 = T02 + m2*x_min_2
f_min_2 = n_2*math.sqrt(gamma*R*T_min_2)/(4*L)

u1 = []
d1 = []
u2 = []
d2 = []

# Finding u from momentum equation
for i in range(len(x)):
    d_temp1 = (P0/(R*T1[i]))
    u1.append((-1*z_min_1[i])/(1j * 2 * math.pi * f_min_1 * d_temp1))   
    d1.append(d_temp1)
    d_temp2 = (P0/(R*T2[i]))
    u2.append((-1*z_min_2[i])/(1j * 2 * math.pi * f_min_2 * d_temp2))   
    d2.append(d_temp2)

# Finding p and u amplitudes 
p_abs_1 =  [abs(x) for x in min_P_1]
p_abs_2 =  [abs(x) for x in min_P_2]
u_abs_1 =  [abs(x) for x in u1]
u_abs_2 =  [abs(x) for x in u2]

plt.plot(x,p_abs_1)
plt.plot(x,p_abs_2)
plt.legend(["T0 = 500K, m = -50, n = 5","T0 = 1100K, m = -200, n = 5"], loc="best")
plt.xlabel("Distance (m)")
plt.ylabel("Acoustic Pressure Amplitude (Pa)")
plt.title("Acoustic Pressure Amplitude vs x")
plt.grid()
plt.show()

plt.plot(x,u_abs_1)
plt.plot(x,u_abs_2)
plt.legend(["T0 = 500K, m = -50, n = 5","T0 = 1100K, m = -200, n = 5"], loc="best")
plt.xlabel("Distance (m)")
plt.ylabel("Acoustic Velocity Amplitude (m/s)")
plt.title("Acoustic Velocity Amplitude vs x")
plt.grid()
plt.show()

# Using subplots so that p and u with different orders of magnitude can be plotted on the same plot
fig, ax1 = plt.subplots()

# Plot the p values
ax1.plot(x, p_abs_2, color='r', label='P')
ax1.set_xlabel('Distance (m)')
ax1.set_ylabel('Acoustic Pressure Amplitude (Pa)', color='r')

# Create a secondary y-axis and plot the u values
ax2 = ax1.twinx()
ax2.plot(x, u_abs_2, color='b', label='U')
ax2.set_ylabel('Acoustic Velocity Amplitude (m/s)', color='b')

# Add legends
lines_1, labels_1 = ax1.get_legend_handles_labels()
lines_2, labels_2 = ax2.get_legend_handles_labels()
ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc='best')

plt.title("Nodes and Antinodes")
plt.grid()
plt.show()

##### PART-C of code for comparing  numerical and analytical pressure amplitudes #####

T1 = 1100
T0 = T1
T2 = 300
m = -200
fr = 157.51
w = 2*math.pi*fr
P1 = 2000
gamma = 1.4
R = 287
a = abs(m)*math.sqrt(gamma*R)/2
L = 4
n = 5
h = 0.01
y0 = [2000,0]

# Constants found from BC conditions
c1 = -math.pi*w*math.sqrt(T1)*P1*yn(1,w*math.sqrt(T1)/a)/(2*a)
c2 = math.pi*w*math.sqrt(T1)*P1*jn(1,w*math.sqrt(T1)/a)/(2*a)

x = np.arange(0,L+h,h)

T = []
s = [] # used to transform wave equation to mean temperature space

# Storing T and s values at each location.
for i in range(len(x)):
    T_temp = (T1 + m*x[i])
    s.append(w*math.sqrt(T_temp)/a)
    T.append(T_temp)

p_nu , z = rk_4(f, x, y0, h, m, T0, L, n, w)

p_an = []
u_an = []
u_nu = []

# finding analytical p and u and numerical u
for i in range(len(x)):
    d_temp = (P0/(R*T[i]))
    p_an.append(c1*jn(0,s[i]) + c2*yn(0,s[i]))
    u_an.append((-m/abs(m)) *1j * (c1*jn(1,s[i]) + c2*yn(1,s[i]))/(d_temp*math.sqrt(gamma * R * T[i])))
    u_nu.append((-1*z[i])/(1j * w * d_temp))

# Finding amplitudes
p_abs_nu =  [abs(x) for x in p_nu] 
u_abs_nu =  [abs(x) for x in u_nu]
p_abs_an =  [abs(x) for x in p_an]
u_abs_an =  [abs(x) for x in u_an]

plt.plot(x,p_abs_an,'o',markersize=2.5)
plt.plot(x,p_abs_nu)
plt.legend(["Analytical","Numerical"], loc="best")
plt.xlabel("Distance (m)")
plt.ylabel("Acoustic Pressure Amplitude (Pa)")
plt.title("Acoustic Pressure Amplitude vs x")
plt.grid()
plt.show()

plt.plot(x,u_abs_an,'o',markersize=2.5)
plt.plot(x,u_abs_nu)
plt.legend(["Analytical","Numerical"], loc="best")
plt.xlabel("Distance (m)")
plt.ylabel("Acoustic Velocity Amplitude (Pa)")
plt.title("Acoustic Velocity Amplitude vs x")
plt.grid()
plt.show()
