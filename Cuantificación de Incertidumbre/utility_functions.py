# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:01:08 2022

@author: Luis Medina
"""


import numpy as np

m=1.0
y0=2.0       #[cm]
y0_dot=0.0   #[cm/s]

to=0.0 #[s]
tf=5.0 #[s]
n=500 #number of points
t=np.linspace(to, tf,n)


def y(samples):
     
    """ 
    This function evaluates the free response for a 1 DOF damped system [3]
    note:tt is the variable name for time. It is not term t since t is also
    the "variable name" used by scipy when importing t Student features.
    """ 
   
    c,k=samples[0]
    
    
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


