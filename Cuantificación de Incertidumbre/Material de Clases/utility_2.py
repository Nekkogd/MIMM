# -*- coding: utf-8 -*-
"""
IDMI310-22
Introducción a la cuantificación de incertidumbres en Ingeniería
Magister en Ingeniería Mecánica y Materiales
Facultad de ciencias de la ingeniería
Universidad Austral de Chile
________________________________________
This is a file that contains utility functions.
That means that here may be defined several functions
which can be used (called from)by examples or  other python scripts.
You can add or modify some of those functions as you require 


@author: Luis Medina
"""

"""
Fecha: 09 de Noviembre, 2022
Archivo de definicion de funciones.

** Archivo auxiliar

Contiene funciones definidas, que luego pueden ser llamadas por otro archivo (ex: Example2_2nn.py)

@author: nekko

"""

import numpy as np

# Importing SymPy. SymPy is a Python library for symbolic mathematics
# Find more information about at https://www.sympy.org/en/index.html
from sympy import *

#Defining model system's variables
m,c,k,tt,y0,y0_dot=var('m,c,k,tt,y0,y0_dot')


def y(m,c,k,tt,y0,y0_dot):
    """ 
    This is the model function, that is y=f(q), where q is (are) the parameter(s).
    The function y can return a symbolic expression, being function of the
    corresponding parameters (in this case m,c,k,t,y0 and y0_dot)
    """ 
    # Here, f(q) corresponds to the model defined at Example 2
    wn=sqrt(k/m)          # natural frequency
    zeda=c/(2*sqrt(k*m))  # damping factor
    wd=wn*sqrt(1-zeda**2) # damped natural frequency

    #System's  free response [2]
    #response's amplitude
    Y=sqrt((y0*wn)**2+y0_dot**2+2*zeda*wn*y0*y0_dot)/wd

    #response's phase
    phi=atan(y0*wd/(y0_dot+zeda*wn*y0))

    res=Y*exp(-zeda*wn*tt)*sin(wd*tt+phi)
    
    return res

#sc_c=simplify(diff(y(m,c,k,tt,y0,y0_dot), c))
#sc_k=simplify(diff(y(m,c,k,tt,y0,y0_dot), k))

def Dev_y_respect_to_c(m,c,k,tt,y0,y0_dot):
    """ 
    This is the function' derivative respect to c. It is required to
    calculate the sensitiviy matrix. The derivative has been calculated
    symbollicaly, using Sympy.
    
    """ 
    res=m*(2*c*k*m*sqrt(k/m)*((-c**2 + 4*k*m)/(k*m))**(3/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)*(c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))*sin(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)))) + k**3*m**2*y0*y0_dot*((-c**2 + 4*k*m)/(k*m))**(5/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)*sin(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)))) - k**2*tt*sqrt(k*m)*((-c**2 + 4*k*m)/(k*m))**(5/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)*(c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))*sin(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)))) - sqrt(k/m)*sqrt(k*m)*(-c**2 + 4*k*m)**2*(c*tt*sqrt(k/m)*(k*m)**(3/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2) - 2*k**2*m**2*y0*(c*m*sqrt(k/m)*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)) + k*y0*(-c**2 + 4*k*m)))*(c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))*cos(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))))/(k**3*m**3))*exp(-c*tt*sqrt(k/m)/(2*sqrt(k*m)))/(sqrt(k*m)*sqrt((c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))/(m*sqrt(k*m)))*(-c**2 + 4*k*m)**3*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2))

def Dev_y_respect_to_k(m,c,k,tt,y0,y0_dot):
    """ 
    This is the function' derivative respect to c. It is required to
    calculate the sensitiviy matrix. The derivative has been calculated
    symbollicaly, using Sympy.
    
    """ 
    res=m**2*(-c**2*(k/m)**(3/2)*sqrt(k*m)*((-c**2 + 4*k*m)/(k*m))**(3/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)*(c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))*sin(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)))) + k**3*m**2*y0**2*(k/m)**(3/2)*((-c**2 + 4*k*m)/(k*m))**(5/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)*sin(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)))) - (k/m)**(3/2)*(k*m)**(3/2)*((-c**2 + 4*k*m)/(k*m))**(5/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)*(c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))*sin(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)))) + sqrt(k*m)*(-c**2 + 4*k*m)**2*(c**2*tt*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2) - 8*m*y0*(k*m)**(3/2)*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)) + tt*(-c**2 + 4*k*m)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2))*(c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))*cos(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))))/(2*m**4))*exp(-c*tt*sqrt(k/m)/(2*sqrt(k*m)))/(k**2*sqrt((c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))/(m*sqrt(k*m)))*(-c**2 + 4*k*m)**3*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2))
    return res


def Sensitivity_matrix(m_val,c_val,k_val,t_vec,y0_val,y0_dot_val):
    """ 
    This function evaluates the sensitivity matrix [3] asociated to y=f(q)
    It depends on how y() function has been defined above.
    """ 
    n=len(t_vec)  # is the time vector dimension 
    #Calculating symbolic derivatives of f(q)  respect to q parameter
   
    
    
    sens_vec=np.ndarray((n,2))
       
    for i in range(n):    
           
        sc_num_c=Dev_y_respect_to_c(m_val,c_val,k_val,t_vec[i],y0_val,y0_dot_val)
        sc_num_k=Dev_y_respect_to_k(m_val,c_val,k_val,t_vec[i],y0_val,y0_dot_val)
        sens_vec[i,0]=sc_num_c
        sens_vec[i,1]=sc_num_k
       
    res=np.matmul(np.transpose(sens_vec),sens_vec) #Sensitivy matrix [2]
    
    return res





  


    
# #System's parameter values  [2] 
# #Note: units are assumed since example in ref. [2] does not specified units
# k_val=20.5 #[N/cm]
# m_val=1.0 # [Kg]
# c_val=1.5 # [Ns/cm] <-- used to generate synthetic data

# #initial condition values [2]
# y0_val=2.0       #[cm]
# y0_dot_val=0.0   #[cm/s]

# #Define a time vector to get  system's response and 
# to=0.0 #[s]
# tf=5.0 #[s]
# n=500 #number of points
# t_vec=np.linspace(to, tf,n)

# sigma_error=.1                 # error's standard deviation
# mu_error=0.0                   # error's mean value   
# var_error=sigma_error**2       # error variance



# sen_mat=Sensitivity_matrix(m_val,c_val,k_val,t_vec,y0_val,y0_dot_val)

# V=np.linalg.inv(sen_mat)*var_error                #damping's variance  

   


   
