# -*- coding: utf-8 -*-
"""
IDMI310-22
Introducción a la cuantificación de incertidumbres en Ingeniería
Magister en Ingeniería Mecánica y Materiales
Facultad de ciencias de la ingeniería
Universidad Austral de Chile
________________________________________
Created on Mon Aug 29 16:28:30 2022

References:
[1]«SymPy». https://www.sympy.org/en/index.html (accedido 25 de abril de 2022).
       
    
[2] R. C. Smith, Uncertainty quantification: theory, implementation, and 
    applications, vol. 12. Siam, 2013., Example 7.15, 
    

[3] Singiresu S. Rao. Mechanical vibrations; SI conversion by Philip Griffin.
Harlow, United Kingdom : Pearson, [2018].  



@author: Luis Medina
"""

import numpy as np
from scipy.stats import norm
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt



#System's parameter values  [2] 
#Note: units are assumed since example in ref. [2] does not specified units
k_0=20.5 #[N/cm]
m_0=1.0 # [Kg]
c_0=1.5 # [Ns/cm] <-- used to generate synthetic data

#initial condition values [2]
y0_val=2.0       #[cm]
y0_dot_val=0.0   #[cm/s]

#Define a time vector to get  system's response and 
to=0.0 #[s]
tf=5.0 #[s]
n=500 #number of points
t_vec=np.linspace(to, tf,n)


# Defining functions
#_____________________________________________________________________________
def y(m,c,k,t,y0,y0_dot):
    """ 
    This function evaluates the free response for a 1 DOF damped system [3]
    """ 
    wn=np.sqrt(k/m)          # natural frequency
    zeda=c/(2*np.sqrt(k*m))  # damping factor
    wd=wn*np.sqrt(1-zeda**2) # damped natural frequency

    #System's  free response [2]
    #response's amplitude
    Y=np.sqrt((y0*wn)**2+y0_dot**2+2*zeda*wn*y0*y0_dot)/wd

    #response's phase
    phi=np.arctan(y0*wd/(y0_dot+zeda*wn*y0))

    res=Y*np.exp(-zeda*wn*t)*np.sin(wd*t+phi)
    
    return res






def L(c,m,k,t,y0,y0_dot,y_obs):
    """ 
    This function defines a least squares estimator as function of parameter c [2].
    Least squares estimator is: argmin [sum(Err^2)], where Err= y_obs(i)-y(i).
    y_obs: system'response measurement
    y: system response model 
    
    """ 
    
    # using list comprenhension as an option to avoid a loop
    y_list=[y(m,c,k,t[i],y0,y0_dot) for i in range(n)]
    y_=np.array(y_list) # to convert the list in a numpy array
    #______________________________________________________
    #Alternatively to the loop inside the list cmprenhension,
    #you can use a "classic loop" (less Pythonic but maybe more readable)
    #y_=np.ndarray((n, ))
    # for i in range(n):
        
    #     y_[i]=y(m_0,c_0,k_0,t_vec[i],y0_val,y0_dot_val)
    # It will produce same result than the option used above
    #______________________________________________________   
   
    Err=y_obs-y_
    # Finally the function to be minimized is
    res=np.matmul(Err,np.transpose(Err)) # L=res
       
    return res
 


    






   
#_____________________________________________________________________________


#1 Define noise to be added to sistem´s response. This is a strategy  to 
# generate synthetic data since there is not measurements available

#Assuming that error is iid and#error ~ N(mu_error, var_error ) [2]

sigma_error=.1                 # error's standard deviation
mu_error=0.0                   # error's mean value   
var_error=sigma_error**2       # error variance
error=sigma_error*np.random.randn(len(t_vec))  # errors



#2.System's response "measurments"  
   
y_vec=np.ndarray((n, ))
 
for i in range(n):
    
    y_vec[i]=y(m_0,c_0,k_0,t_vec[i],y0_val,y0_dot_val)
    
   
y_obs=y_vec+error  # "measurements"  



#Find c based on least squares estimator (defined as L function)

m=m_0
k=k_0
t=t_vec
y0=y0_val
y0_dot=y0_dot_val
obs=y_obs
#Define a interval to search the c value that minimizes the Least square estimator
c_min=0.0
c_max=4.0
res_opt=minimize_scalar(L,args=(m,k,t,y0,y0_dot,y_obs), bounds=(c_min, c_max), method='bounded')
c_opt=res_opt.x
print('El c optimo es:', c_opt)



   
#Plotting results
fig, axs=plt.subplots(2,1)
axs[0].plot(t,y_vec, linewidth=1)
axs[0].scatter(t,y_obs,s=0.1)

axs[0].set_xlabel('t [s]')
axs[0].set_ylabel('y [cm]')
axs[0].grid(True)

axs[1].scatter(t,error,s=1)
axs[1].plot(t,-2*sigma_error*np.ones(len(t)),"b-",linewidth=1)
axs[1].plot(t,2*sigma_error*np.ones(len(t)),"b-",linewidth=1)
axs[1].set_xlabel('t [s]')
axs[1].set_ylabel('error [cm]')
axs[1].grid(True)
fig.tight_layout()
plt.show()