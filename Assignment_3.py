#Importing necessary modules
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from scipy.signal import argrelextrema
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D
from matplotlib.cm import ScalarMappable
from matplotlib.gridspec import GridSpec
from matplotlib.colors import ListedColormap



#Function to return eta and eta_dot to use in rk4 for t < t_d
def f1(y, j, c1, c2, K_i, t_d, M,x_f):
    gamma = 1.4
    
    a, b = y  
    k_j = j*math.pi
    w_j = k_j
    #z_j = 0   # for no damping case
    z_j = (c1*w_j/math.pi + c2*math.sqrt(math.pi/w_j))/(2*math.pi)     # damping

    db_dt = -2*z_j*w_j*b - (k_j**2)*a 
    
    return [b, db_dt]

#Function to return eta and eta_dot to use in rk4 for t > t_d
def f2(y, j, c1, c2, K_i, t_d, M, x_f, u_f):
    gamma = 1.4
    
    a, b = y  
    k_j = j*math.pi
    w_j = k_j
    #z_j = 0  # for no damping case
    z_j = (c1*w_j/math.pi + c2*math.sqrt(math.pi/w_j))/(2*math.pi) # damping

    db_dt = -2*z_j*w_j*b - (k_j**2)*a - j*math.pi*K_i*(math.sqrt(abs(1/3 + u_f)) - math.sqrt(1/3))*math.sin(j*math.pi*x_f)

    return [b, db_dt]

# Function to find eta and eta for t < t_d
def rk_4_1(f1, y0, h, j, c1, c2, K_i, t_d, M, x_f):
    a, b = [y0[0]], [y0[1]]
   
    for i in range(int(t_d/h)):
        a_i, b_i  = a[-1], b[-1] 
        k1 =  h*np.array(f1([a_i, b_i], j, c1, c2, K_i, t_d, M, x_f))
        k2 =  h*np.array(f1([a_i + (k1[0] / 2), b_i + (k1[1] / 2)], j, c1, c2, K_i, t_d, M, x_f))
        k3 =  h*np.array(f1([a_i + (k2[0] / 2), b_i + (k2[1] / 2)], j, c1, c2, K_i, t_d, M, x_f))
        k4 =  h*np.array(f1([a_i + k3[0], b_i + k3[1]], j, c1, c2, K_i, t_d, M, x_f))
        a_i_1 = a_i + (1 / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
        b_i_1 = b_i + (1 / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
        a.append(a_i_1)
        b.append(b_i_1)

    return a, b

# Function to find eta and eta for t > t_d
def rk_4_2(f1, f2, t, y0, h, c1, c2, K_i, t_d, M, x_f, N):
    
    eta = []
    eta_dot = []
    
    l = int(t_d/h) + 1

    for i in range(N):
        eta_j, eta_dot_j = rk_4_1(f1, y0[i], h, i+1, c1, c2, K_i, t_d, M, x_f)
        eta.append(eta_j)
        eta_dot.append(eta_dot_j)

    u1 = []                   
    for k in range(l):
        u_temp = 0
        for i in range(N):
            u_temp += eta[i][k]*math.cos((i+1)*math.pi*x_f)
        u1.append(u_temp)

   
    for i in range(len(t) - l):    # marching in time
        u_temp = 0
        for k in range(N):    # amrching in modes

            a_1 = eta[k]  
            b_1 = eta_dot[k]
            a_i, b_i  = a_1[-1], b_1[-1]

            k1 =  h*np.array(f2([a_i, b_i], k+1, c1, c2, K_i, t_d, M, x_f, u1[i+1]))
            k2 =  h*np.array(f2([a_i + (k1[0] / 2), b_i + (k1[1] / 2)], k+1, c1, c2, K_i, t_d, M, x_f, u1[i+1]))
            k3 =  h*np.array(f2([a_i + (k2[0] / 2), b_i + (k2[1] / 2)], k+1, c1, c2, K_i, t_d, M, x_f, u1[i+1]))
            k4 =  h*np.array(f2([a_i + k3[0], b_i + k3[1]], k+1, c1, c2, K_i, t_d, M, x_f, u1[i+1]))
            a_i_1 = a_i + (1 / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0])
            b_i_1 = b_i + (1 / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
            u_temp += a_i_1 * math.cos((k+1)*math.pi*x_f)               # summing up eta_jcos(j*pi*x_f) at a particular time for all modes
            a_1.append(a_i_1)
            b_1.append(b_i_1)
            eta[k] = a_1
            eta_dot[k] = b_1
    
        u1.append(u_temp)       # appending the summed up u' for each time
    
    return eta, eta_dot, u1

#Initial values found from references. Change K, t_d, T and N and x_f according to requirements
gamma = 1.4

T = 80
h = 0.1
c_0 = 399.6
L_w = 3.6
u_bar = 0.5
M = u_bar/c_0
rho_bar = 1.205
T_w = 1000
C_v = 719
lamda = 0.0328
T_bar = 300
d_w = 0.0005
S = 0.0018
p_bar = 101325
c1 = 0.1
c2 = 0.06

x_f = 0.29  # flame location
t = np.arange(0,T+h,h)
K = 1  # Heater power
#K =  (gamma - 1)*2*L_w*(T_w - T_bar)*math.sqrt(math.pi*lamda*C_v*rho_bar*u_bar*d_w/2)/(S*c_0*p_bar*math.sqrt(3))  # formula for K, gives decayng solution

t_d = 0.5   # time delay
N = 3
n = []

for i in range(1, N+1):
    n.append(i)

y0 = [[0.18, 0]] + [[0, 0]] * (N - 1)   # conditions of eta and eta_dot at t = 0

################## PART OF CODE TO REPRODUCE FIG 4, 5 AND 6 PF REFERENCE 2 ##################
eta, eta_dot, u = rk_4_2(f1, f2, t, y0, h, c1, c2, K, t_d, M, x_f, N)
p = []

for k in range(len(t)):  # 0,1,2,
        p_temp = 0
        for i in range(N): 
            p_temp += (gamma*M*eta_dot[i][k]*math.sin(n[i]*math.pi*x_f))/(n[i]*math.pi)
        
        p.append(p_temp)

E = []

for i in range(len(u)):
    E.append((0.5*(p[i]**2) + 0.5*((gamma*M*u[i])**2))/((gamma*M)**2))

#print(eta[1])
#print(eta[2])

#print(len(u))

# Find the indices of local maxima (peaks)
peaks_idx = argrelextrema(np.array(E), np.greater, order=7)[0]

# Filter peaks to keep only the highest ones
min_prominence = 0  # Adjust as needed
prominences = np.array([E[idx] for idx in peaks_idx])
highest_peaks_idx = peaks_idx[prominences > min_prominence]

plt.plot(t, u)
plt.title("Non-dimensional velocity vs time")
plt.xlabel("t")
plt.ylabel("u'")
#plt.ylim([-0.2,0.2])
plt.grid()
plt.show()

plt.plot(t,E, label='Energy')
plt.plot(np.array(t)[highest_peaks_idx], np.array(E)[highest_peaks_idx], color='red', markersize=10, label='Peaks')
plt.title("Non-Linear Evolution of Acoustic Energy")
plt.xlabel("t")
plt.ylabel("$\\frac{Energy}{\\gamma M^{2}}$")
plt.legend()
plt.grid()
plt.show()

plt.plot(t,eta[0])
plt.title("Acoustic velocity projected onto first Galerkin Mode")
plt.xlabel("t")
plt.ylabel("$\\eta_1}$")
plt.grid()
plt.show()

plt.plot(t,eta[1])
plt.title("Acoustic velocity projected onto second Galerkin Mode")
plt.xlabel("t")
plt.ylabel("$\\eta_2}$")
plt.grid()
plt.show()

plt.plot(t,eta[2])
plt.title("Acoustic velocity projected onto third Galerkin Mode")
plt.xlabel("t")
plt.ylabel("$\\eta_3}$")
plt.ylim([-0.2,0.2])
plt.grid()
plt.show()

################## PART OF CODE TO REPRODUCE FIG 6A (2D BIFURCATION PLOT) OF  REFERENCE 1 ##################
T = 100
x_f = 0.3
t_d = 0.5
h = 0.1
t = np.arange(0,T+h,h)
K = np.linspace(0.05, 2, 30)
K_rev = np.linspace(2,0.05, 30)

U1 = []
U2 = []
y0 = [[0.18, 0]] + [[0, 0]] * (N - 1)

# Increasing K
for i in range(len(K)):
    a_2, b_2, u_f = rk_4_2(f1, f2, t, y0, h, c1, c2, K[i], t_d, M, x_f, N)
    second_half_index = len(u_f) // 2
    U1.append(np.sqrt(np.mean(np.square(u_f[second_half_index:])))) # Using second half of u' values to calculate Urms

    for x in range(N):
        y0[x] = [a_2[x][len(t)-1] , b_2[x][len(t)-1]]   # Using the eta_j(T) and eta_dot_j(T) of ith K as the initial condition for (i+1)th K.

# Decreasing K
for i in range(len(K_rev)):
    a_2_rev, b_2_rev, u_f_rev = rk_4_2(f1, f2, t, y0 , h, c1, c2, K_rev[i], t_d, M, x_f, N)
    second_half_index = len(u_f) // 2
    U2.append(np.sqrt(np.mean(np.square(u_f_rev[second_half_index:]))))

    for x in range(N):    
        y0[x] = [a_2[x][len(t)-1] , b_2[x][len(t)-1]]


plt.plot(K,U1,'o',label = "Increasing K")
plt.plot(K_rev,U2,'*', label = "Decreasing K")
plt.xlabel('K')
plt.ylabel('$|U_{1}|$')
plt.title('Bifurcation plot with K as parameter')
plt.legend()
plt.grid()
plt.show()

################## PART OF CODE TO REPRODUCE FIG 6B (3D BIFURCATION PLOT) OF  REFERENCE 1 ##################

T = 80
x_f = 0.29
N = 1
h = 0.01
t = np.arange(0,T+h,h)
K = np.linspace(0.4, 2, 50)
K_rev = np.linspace(2, 0.4, 50)
t_d = np.linspace(0.2, 0.6, 50)

U1 = np.zeros((len(K), len(t_d)))
U2 = np.zeros((len(K), len(t_d)))

for k in range(len(t_d)):
    y0 = [[0.2, 0]] + [[0, 0]] * (N - 1)
    for i in range(len(K)):
   
        a_2, b_2, u_f = rk_4_2(f1, f2, t, y0, h, c1, c2, K[i], t_d[k], M, x_f, N)

        second_half_index = len(u_f) // 2
        U1[i, k] = np.sqrt(np.mean(np.array(u_f[second_half_index:]) ** 2))
       
        for x in range(N):
            y0[x] = [a_2[x][len(t)-1] , b_2[x][len(t)-1]]

    for i in range(len(K)):
        a_2_rev, b_2_rev, u_f_rev = rk_4_2(f1, f2, t, y0 , h, c1, c2, K_rev[i], t_d[k], M, x_f, N)
        second_half_index = len(u_f) // 2
        U2[i,k] = (np.sqrt(np.mean(np.square(u_f_rev[second_half_index:]))))

        for x in range(N):    
            y0[x] = [a_2[x][len(t)-1] , b_2[x][len(t)-1]]

# Creating mesh
K_mesh, t_d_mesh = np.meshgrid(K, t_d)
K_mesh_rev, t_d_mesh_rev = np.meshgrid(K_rev, t_d)

# Plotting
fig = plt.figure(figsize=(16, 10))
cmap_single_blue = ListedColormap(['blue'])
# 3D plot
ax1 = fig.add_subplot(111, projection='3d')
surf1 = ax1.plot_surface(K_mesh, t_d_mesh, U1.T, cmap = cmap_single_blue)
surf2 = ax1.plot_surface(K_mesh_rev, t_d_mesh_rev, U2.T, cmap='viridis')
ax1.set_title('$|U_{1}|$ vs $\\tau$ vs K')
ax1.set_xlabel('K')
ax1.set_ylabel('$\\tau$')
ax1.set_zlabel('$|U_{1}|$')

# Colorbars
cb1 = plt.colorbar(surf1, shrink=0.3, aspect= 20, pad=0.0)
cb1.set_label('Increasing K')
cb2 = plt.colorbar(surf2, shrink=0.3, aspect=20, pad=0.001)
cb2.set_label('Decreasing K')

plt.show()
