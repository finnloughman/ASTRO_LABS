# -*- coding: utf-8 -*-

"""
Created on Fri Sep 30 18:13:58 2022

@author: finnl
"""

import numpy as np
import matplotlib.pyplot as plt

print('Finn Loughman - 20332121')


#A
z = ([[0 for i in range(2)] for x in range(10)])


#B
def func(z,a,b):
    return sum(z)**a, sum(z)**b

rand = np.random.randint(0,2,10)
r_array = np.array(rand)

print(func(r_array,2,3))

x = func(r_array,2,3)[0]
y = func(r_array,2,3)[1]
print (x,y)

#C
a = np.random.normal(10,3,500)

plt.hist(a,20,edgecolor='b')
plt.xlabel('Number')
plt.ylabel('Frequency')
plt.title('Gaussian Distribution of 500 random numbers to 3 standard deviations')
plt.show()

#D
c = np.array([-10,10])
b = np.random.choice(c,20)

for x in range(3):
    print (b)
  
#%%
#E
y1 = []

d = np.arange(-50,51,1)

for i in range(0,50):
        count =0
        while count<100:
            count += 1
            y1.append(np.random.randint(-50,50))
                

plt.hist(y1,30,alpha=0.7,edgecolor='b')
plt.xlabel('Number')
plt.ylabel('Frequency')
plt.title('Distribution of 5000 random numbers between -50 and 50')

"""
Created on Thu Nov  3 23:13:15 2022

@author: finnl
"""
import numpy as np
import matplotlib.pyplot as plt


print('Finn Loughman - 20332121')
#AAAAAAAAAAAAAAAAAAAAAAAA

#Defining Variables
f = 10
k = 1 
M = 100
m = 1 
T_0 = 1
sigma_par = np.sqrt(k*T_0/m) 
sigma_pis = np.sqrt(k*T_0/M)
mu = 0



#Positions and Velocities  
x_pis = 100     #(ideal gas law gives eq position N/F = 100)
x_par = np.random.uniform (low=-x_pis , high = x_pis, size = 1000)
v_par = np.random.normal(mu, sigma_par, 1000)
v_pis = np.random.normal(mu, sigma_pis, 1)


#Waiting Time tau
def Tau(M, f, v_pis, v_par, x_pis, x_par):

    if v_par > 0:
        y =(M/f) * (v_pis - v_par + np.sqrt((v_pis-v_par)**2-2*(f/M)*(x_par-x_pis)))  #Right moving particles
        return y

    else:
        y = (M/f) * (v_pis + v_par + np.sqrt((v_pis+v_par)**2+2*(f/M)*(x_par+x_pis))) #Left moving particles
        return y
    
#Initial Setup
Tau = np.vectorize(Tau)
times0 = []
times1 = []
pos_pist = []
t_step = []
v_pist_list = []

H_list = []




for n in range(0, 20000):
    t_step = Tau(M, f, v_pis, v_par, x_pis, x_par)
    time_step = min(t_step)
    t_step = t_step.tolist()
    v_index = t_step.index(time_step)
    x_par += v_par * time_step
    x_pis = abs(x_par[v_index])
    v_pis -= (f/M) * time_step
    vs = v_par[v_index]
    h = sum(0.5*m*v_par**2) +0.5*M*(v_pis**2) +f*x_pis
     
    if vs > 0: #right moving
        v_pis_new = ((2*m*vs+(M-m)*v_pis)/(M+m))
        v_par[v_index] = ((2*M*v_pis+(m-M)*vs)/(M+m)) 
    else:  #left moving
        v_pis_new = ((-2*m*vs+(M-m)*v_pis)/(M+m))
        v_par[v_index] = ((-2*M*v_pis+(m-M)*vs)/(M+m))
        
    v_pis = v_pis_new
    times0.append(time_step)
    z = sum(times0)
    pos_pist.append(x_pis)
    times1.append(z)
    v_pist_list.append(v_pis)
    H_list.append(h)
    n += 1


plt.plot(times1,pos_pist,label='Position(t)')
plt.xlabel("Time")
plt.ylabel("Piston Position")
plt.title("Piston Position over time")
plt.axhline(y=np.average(pos_pist), color ='r', linestyle ='dashed',label='Expected Position')
plt.legend(loc=1)
plt.show()

plt.plot(times1, v_pist_list, label='Velocity(t)', linewidth=0.5)
plt.axhline(y = np.average(v_par), color = 'r', linestyle = 'dashed',label='Average Velocity')
plt.xlabel("Time")
plt.ylabel("Piston Velocity")
plt.title("Piston Velocity over Time")
plt.legend(loc=1)
plt.show()

plt.plot(times1, H_list)
plt.xlabel("Time")
plt.ylabel("Enthalpy")
plt.title("Enthalpy over time")
plt.show()

F_list=[0.1,0.25,1,5,10,50,100]

weighted_av= np.average(pos_pist[4500:])
print(weighted_av)

#%%BBBBBBBBBBBBBBBBBBBB
F_list=[0.1,0.25,1,5,10,50,100]
average_pos =[]
for f in F_list:
    k = 1 
    M = 100
    m = 1 
    T_0 = 1
    sigma_par = np.sqrt(k*T_0/m) 
    sigma_pis = np.sqrt(k*T_0/M)
    mu = 0
    
    
    
    #Positions and Velocities  
    x_pis = 100     #(ideal gas law gives eq position N/F = 100)
    x_par = np.random.uniform (low=-x_pis , high = x_pis, size = 1000)
    v_par = np.random.normal(mu, sigma_par, 1000)
    v_pis = np.random.normal(mu, sigma_pis, 1)
    
    
    #Waiting Time tau
    def Tau(M, f, v_pis, v_par, x_pis, x_par):
    
        if v_par > 0:
            y =(M/f) * (v_pis - v_par + np.sqrt((v_pis-v_par)**2-2*(f/M)*(x_par-x_pis)))  #Right moving particles
            return y
    
        else:
            y = (M/f) * (v_pis + v_par + np.sqrt((v_pis+v_par)**2+2*(f/M)*(x_par+x_pis))) #Left moving particles
            return y
        
    #Initial Setup
    Tau = np.vectorize(Tau)
    times0 = []
    times1 = []
    pos_pist = []
    t_step = []
    v_pist_list = []
    
    H_list = []
    
    
    
    
    for n in range(0, 20000):
        t_step = Tau(M, f, v_pis, v_par, x_pis, x_par)
        time_step = min(t_step)
        t_step = t_step.tolist()
        v_index = t_step.index(time_step)
        x_par += v_par * time_step
        x_pis = abs(x_par[v_index])
        v_pis -= (f/M) * time_step
        vs = v_par[v_index]
        h = sum(0.5*m*v_par**2) +0.5*M*(v_pis**2) +f*x_pis
         
        if vs > 0: #right moving
            v_pis_new = ((2*m*vs+(M-m)*v_pis)/(M+m))
            v_par[v_index] = ((2*M*v_pis+(m-M)*vs)/(M+m)) 
        else:  #left moving
            v_pis_new = ((-2*m*vs+(M-m)*v_pis)/(M+m))
            v_par[v_index] = ((-2*M*v_pis+(m-M)*vs)/(M+m))
            
        v_pis = v_pis_new
        times0.append(time_step)
        z = sum(times0)
        pos_pist.append(x_pis)
        times1.append(z)
        v_pist_list.append(v_pis)
        H_list.append(h)
        n += 1
    
    
    weighted_av= np.average(pos_pist[4500:],weights=times1)
    average_pos.append(weighted_av)
    
ideal_gas = []
for f in F_list:
    x =1000/f
    ideal_gas.append(x)
    
plt.scatter(F_list,average_pos,label= 'System data points',c='r')
plt.xscale('log')
plt.yscale('log')
plt.loglog(F_list,ideal_gas,label='Ideal gas curve')
plt.title('Weighted average of pistion position as a function of force')
plt.xlabel('Forces')
plt.ylabel('Average piston position')
plt.legend(loc=1)
plt.show()

#%% CCCCCCCCCCCCCCCCC
import scipy.stats
from scipy.stats import maxwell

print('Finn Loughman - 20332121')

f = 10
k = 1 #.38e-23
M = 100 #*1.66e-27
m = 1 #*1.66e-27
T_0 = 1
sigma_par = np.sqrt(k*T_0/m) 
sigma_pis = np.sqrt(k*T_0/M)
mu = 0
R= 8.314
v = v_par
#def maxwell(v):
   # return 4*np.pi((M/(2*np.pi*R*T_0))**(3/2))*(v**2)*np.exp((-M*v**2)/2*R*T_0)
#Position of piston   (ideal gas law gives eq position N/F = 100)
x_pis = 200
#Position of particle
x_par = np.random.uniform(low=-x_pis , high = x_pis, size = 1000)
#Velocity of particle
# = np.random.normal(mu,sigma,1)[0]
#rand_array = np.array([v_part_random])
v_rand =[-5.2,5.2]
v_par = (np.random.choice(v_rand,1000))
#velocity of piston
v_pis = np.random.normal(mu, sigma_pis, 1)

def t(M, f, v_pis, v_par, x_pis, x_par):
    
#Right moving particles
    if v_par > 0:
        y =(M/f) * (v_pis - v_par + np.sqrt((v_pis-v_par)**2-2*(f/M)*(x_par-x_pis)))
        return y
#Left moving particles
    else:
        y = (M/f) * (v_pis + v_par + np.sqrt((v_pis+v_par)**2+2*(f/M)*(x_par+x_pis)))
        return y
    
#Initial Setup
t = np.vectorize(t)
times0 = []
times1 = []
pos_pist = []
t_step = []
v_pist_list = []
new_velocities = []
v_part_list = []
H_list = []



for n in range(0, 20000):
    t_step = t(M, f, v_pis, v_par, x_pis, x_par)
    time_step = min(t_step)
    t_step = t_step.tolist()
    v_index = t_step.index(time_step)
    x_par += v_par * time_step
    x_pis = abs(x_par[v_index])
    v_pis -= (f/M) * time_step
    vs = v_par[v_index]
     
    if vs > 0: #right moving
        v_pis_new = ((2*m*vs+(M-m)*v_pis)/(M+m))
        v_par[v_index] = ((2*M*v_pis+(m-M)*vs)/(M+m)) 
    else:  #left moving
        v_pis_new = ((-2*m*vs+(M-m)*v_pis)/(M+m))
        v_par[v_index] = ((-2*M*v_pis+(m-M)*vs)/(M+m))
        
    v_pis = v_pis_new
    times0.append(time_step)
    z = sum(times0)
    pos_pist.append(x_pis)
    times1.append(z)
    v_pist_list.append(v_pis)
    new_velocities.append(v_index)
    v_part_list.append(vs)
   # H_list.append(h)
    n += 1

x = np.linspace(-20,20,1000)

plt.hist(v_part_list[:500],bins=30, color = 'b',edgecolor='k')
curve = maxwell.fit(v_part_list[0:500])
plt.plot(x,maxwell.pdf(x,*curve)*3000)
plt.title('0-500 iterations')
plt.xlabel('Particle Velocities')
plt.ylabel('Frequency')
plt.show()
plt.hist(v_part_list[500:2500],bins=40, color = 'r',edgecolor='k')
plt.plot(x,maxwell.pdf(x,*curve)*6000)
plt.title('500-2500 iterations')
plt.xlabel('Particle Velocities')
plt.ylabel('Frequency')
plt.show()
plt.hist(v_part_list[2500:7500],bins=40, color = 'g',edgecolor='k')
plt.plot(x,maxwell.pdf(x,*curve)*5000)
plt.title('2500-7500 iterations')
plt.xlabel('Particle Velocities')
plt.ylabel('Frequency')
plt.show()
plt.hist(v_part_list[7500:12500],bins=40, color= 'orange',edgecolor='k')
plt.plot(x,maxwell.pdf(x,*curve)*4000)
plt.title('7500-12500 iterations')
plt.xlabel('Particle Velocities')
plt.ylabel('Frequency')
plt.show()
plt.hist(v_part_list[12500:15000],bins=40, color= 'purple',edgecolor='k')
plt.plot(x,maxwell.pdf(x,*curve)*2000)
plt.title('12500-15000 iterations')
plt.xlabel('Particle Velocities')
plt.ylabel('Frequency')
plt.show()
plt.hist(v_part_list[15000:17500],bins=40, color= 'c',edgecolor='k')
plt.plot(x,maxwell.pdf(x,*curve)*2000)
plt.title('15000-17500 iterations')
plt.xlabel('Particle Velocities')
plt.ylabel('Frequency')
plt.show()
plt.hist(v_part_list[17500:],bins=40, color='yellow',edgecolor='k')
plt.plot(x,maxwell.pdf(x,*curve)*2000)
plt.title('17500-20000 iterations')
plt.xlabel('Particle Velocities')
plt.ylabel('Frequency')
plt.show()


#%% 3Dplot Histograms
import mpl_toolkits.mplot3d
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
nbins = 50


for c, z , ys in zip(['r', 'b', 'g', 'c', 'm', 'y', 'k'],[250,1000,2250,3750,5250,6750,8250], [v_part_list[:500],v_part_list[500:1500],v_part_list[1500:3000],v_part_list[3000:4500],v_part_list[4500:6000],v_part_list[6000:7500],v_part_list[7500:9000]]):
    
    hist, bins = np.histogram(ys, bins=nbins)
    xs = (bins[:-1] + bins[1:])/2

    ax.bar(xs, hist, zs=z, zdir='y', color=c, ec=c, alpha=0.8)

ax.set_xlabel('Particle Velocities')
ax.set_ylabel('Interation no.')
ax.set_zlabel('Frequency')

plt.show()
