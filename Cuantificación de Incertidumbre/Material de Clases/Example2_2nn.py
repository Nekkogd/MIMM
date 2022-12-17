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
"""
Fecha: 09 de Noviembre, 2022
inslatalamos SymPy en nuestro enviorment antes de empezar

@author: nekko

"""

#Basic libraries
import numpy as np
import matplotlib.pyplot as plt
#For statistics purposes
from scipy.stats import norm
from scipy.stats import t 
from scipy.optimize import minimize_scalar

# an library created for this example
from utility_2 import *

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
def y(m,c,k,tt,y0,y0_dot):
    """ 
    This function evaluates the free response for a 1 DOF damped system [3]
    note:tt is the variable name for time. It is not term t since t is also
    the "variable name" used by scipy when importing t Student features.
    """ 
    wn=np.sqrt(k/m)          # natural frequency
    zeda=c/(2*np.sqrt(k*m))  # damping factor
    wd=wn*np.sqrt(1-zeda**2) # damped natural frequency

    #System's  free response [2]
    #response's amplitude
    Y=np.sqrt((y0*wn)**2+y0_dot**2+2*zeda*wn*y0*y0_dot)/wd

    #response's phase
    phi=np.arctan(y0*wd/(y0_dot+zeda*wn*y0))

    res=Y*np.exp(-zeda*wn*tt)*np.sin(wd*tt+phi)
    
    return res





def L(c,m,k,tt,y0,y0_dot,y_obs):
    """ 
    This cost function defines a least squares estimator as function of parameter c [2].
    Least squares estimator is: argmin [sum(Err^2)], where Err= y_obs(i)-y(i).
    y_obs: system'response measurement
    y: system response model 
    
    """ 
    
    # using list comprenhension as an option to avoid a loop
    y_list=[y(m,c,k,tt[i],y0,y0_dot) for i in range(n)]
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
error_var=sigma_error**2       # error variance
error=sigma_error*np.random.randn(n)  # errors



#2.System's response "measurments"  
   
y_vec=np.ndarray((n, ))
 
for i in range(n):
    
    y_vec[i]=y(m_0,c_0,k_0,t_vec[i],y0_val,y0_dot_val)
    
   
y_obs=y_vec+error  # "measurements"  




#3 Define a interval to search the c value that minimizes the Least square estimator
c_min=0.0
c_max=3.0
res_opt=minimize_scalar(L,args=(m_0,k_0,t_vec,y0_val,y0_dot_val,y_obs), bounds=(c_min, c_max), method='bounded')
c_opt=res_opt.x


#4 Calculating the confidence interval for the damping coefficient estimate
# Calculating the Sensitivity matrix
p=1   #p is the number of parameters to be estimated. In this example p=1
      #since it is the damping coefficient the only parameter to be estimated

#Using a function defined at utility_2 library to get the sensitivity matrix
Sens_matrix=Sensitivity_matrix(m_0,c_opt,k_0,t_vec,y0_val,y0_dot_val)

# Sensitivity matrix is a pxp matrix. Hence, for this example is 1x1 (a scalar)

#Computing an estimation of error variance 
      
est_error_var=np.var(error,ddof=p) # Notice that this value should be close to sigma_error**2
                                   # that is the error variance used to generate
                                   # the y_obs data
                              
                              
est_error_std=np.sqrt(est_error_var) # estimated error standard deviation
                                     # it must be close to sigma_error
#Estimating the damping coeff 's variance 
V=np.linalg.inv(Sens_matrix)*est_error_var                

#Confidence interval for the damping coefficient estimate

alpha=0.95 # confidence level

SE=est_error_std*np.sqrt(Sens_matrix) #standard error

#Finally, with a "alpha" level confidence the min and max values for the
# estimated damping value is

c_min,c_max=t.interval(alpha, n-p, loc=c_opt, scale=SE)

#___________________________________________________________________________   
#Plotting results
fig, axs=plt.subplots(2,1)
axs[0].plot(t_vec,y_vec, linewidth=1)
axs[0].scatter(t_vec,y_obs,s=0.1)

axs[0].set_xlabel('t [s]')
axs[0].set_ylabel('y [cm]')
axs[0].grid(True)

axs[1].scatter(t_vec,error,s=1)
axs[1].plot(t_vec,-2*sigma_error*np.ones(len(t_vec)),"b-",linewidth=1)
axs[1].plot(t_vec,2*sigma_error*np.ones(len(t_vec)),"b-",linewidth=1)
axs[1].set_xlabel('t [s]')
axs[1].set_ylabel('error [cm]')
axs[1].grid(True)
fig.tight_layout()
plt.show()