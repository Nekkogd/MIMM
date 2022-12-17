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


 #Calculating symbolic derivatives of f(q)  respect to q parameters 

def Dev_y_respect_to_c(m,c,k,tt,y0,y0_dot):
    """ 
    This is the function' derivative respect to c. It is required to
    calculate the sensitiviy matrix. The derivative has been calculated
    symbollicaly, using Sympy.
    
    """ 
    res=m*(2*c*k*m*sqrt(k/m)*((-c**2 + 4*k*m)/(k*m))**(3/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)*(c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))*sin(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)))) + k**3*m**2*y0*y0_dot*((-c**2 + 4*k*m)/(k*m))**(5/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)*sin(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)))) - k**2*tt*sqrt(k*m)*((-c**2 + 4*k*m)/(k*m))**(5/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)*(c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))*sin(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)))) - sqrt(k/m)*sqrt(k*m)*(-c**2 + 4*k*m)**2*(c*tt*sqrt(k/m)*(k*m)**(3/2)*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2) - 2*k**2*m**2*y0*(c*m*sqrt(k/m)*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m)) + k*y0*(-c**2 + 4*k*m)))*(c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))*cos(tt*sqrt(k/m)*sqrt((-c**2 + 4*k*m)/(k*m))/2 + atan(y0*sqrt(k/m)*sqrt(k*m)*sqrt((-c**2 + 4*k*m)/(k*m))/(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))))/(k**3*m**3))*exp(-c*tt*sqrt(k/m)/(2*sqrt(k*m)))/(sqrt(k*m)*sqrt((c*m*y0*y0_dot*sqrt(k/m) + k*y0**2*sqrt(k*m) + m*y0_dot**2*sqrt(k*m))/(m*sqrt(k*m)))*(-c**2 + 4*k*m)**3*(k*y0**2*(c**2 - 4*k*m) - m*(c*y0*sqrt(k/m) + 2*y0_dot*sqrt(k*m))**2)) 
    return res
def Dev_y_respect_to_k(m,c,k,tt,y0,y0_dot): 
    
   res=-c**2*sqrt(c*y0*y0_dot*sqrt(k/m)/sqrt(k*m) + k*y0**2/m + y0_dot**2)*exp(-c*tt*sqrt(k/m)/(2*sqrt(k*m)))*sin(tt*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1) + atan(y0*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1)/(c*y0*sqrt(k/m)/(2*sqrt(k*m)) + y0_dot)))/(8*k**2*m*sqrt(k/m)*(-c**2/(4*k*m) + 1)**(3/2)) + sqrt(c*y0*y0_dot*sqrt(k/m)/sqrt(k*m) + k*y0**2/m + y0_dot**2)*(c**2*tt*sqrt(k/m)/(8*k**2*m*sqrt(-c**2/(4*k*m) + 1)) + (c**2*y0*sqrt(k/m)/(8*k**2*m*sqrt(-c**2/(4*k*m) + 1)*(c*y0*sqrt(k/m)/(2*sqrt(k*m)) + y0_dot)) + y0*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1)/(2*k*(c*y0*sqrt(k/m)/(2*sqrt(k*m)) + y0_dot)))/(k*y0**2*(-c**2/(4*k*m) + 1)/(m*(c*y0*sqrt(k/m)/(2*sqrt(k*m)) + y0_dot)**2) + 1) + tt*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1)/(2*k))*exp(-c*tt*sqrt(k/m)/(2*sqrt(k*m)))*cos(tt*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1) + atan(y0*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1)/(c*y0*sqrt(k/m)/(2*sqrt(k*m)) + y0_dot)))/(sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1)) + y0**2*exp(-c*tt*sqrt(k/m)/(2*sqrt(k*m)))*sin(tt*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1) + atan(y0*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1)/(c*y0*sqrt(k/m)/(2*sqrt(k*m)) + y0_dot)))/(2*m*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1)*sqrt(c*y0*y0_dot*sqrt(k/m)/sqrt(k*m) + k*y0**2/m + y0_dot**2)) - sqrt(c*y0*y0_dot*sqrt(k/m)/sqrt(k*m) + k*y0**2/m + y0_dot**2)*exp(-c*tt*sqrt(k/m)/(2*sqrt(k*m)))*sin(tt*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1) + atan(y0*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1)/(c*y0*sqrt(k/m)/(2*sqrt(k*m)) + y0_dot)))/(2*k*sqrt(k/m)*sqrt(-c**2/(4*k*m) + 1)) 
   return res

#Sensitivity matrix (when the q parameter is c)
def Sensitivity_matrix(m_val,c_val,k_val,t_vec,y0_val,y0_dot_val):
    """ 
    This function evaluates the sensitivity matrix [3] asociated to y=f(q)
    It depends on how y() function has been defined above.
    """ 
    n=len(t_vec)  # is the time vector dimension 
    #Calculating symbolic derivatives of f(q)  respect to q parameter
   
    
    
    sens_vec=np.ndarray((n,1))
       
    for i in range(n):    
           
        sc_num=Dev_y_respect_to_c(m_val,c_val,k_val,t_vec[i],y0_val,y0_dot_val)
        
        sens_vec[i,0]=sc_num
       
       
    res=np.matmul(np.transpose(sens_vec),sens_vec) #Sensitivy matrix [2]
    
    return res

#Sensitivity matrix (when the q parameters are c and k)
def Sensitivity_matrix2(m_val,c_val,k_val,t_vec,y0_val,y0_dot_val):
    """ 
    This function evaluates the sensitivity matrix [3] asociated to y=f(q)
    It depends on how y() function has been defined above. In this case is a 2D 
    matrix, since comprises the sensitivities respect to c and k coefficients
    """ 
    n=len(t_vec)  # is the time vector dimension 
   
   
    
    
    sens_vec=np.ndarray((n,2))
       
    for i in range(n):    
        #Evaluating f derivative respect to c   
        sc_num=Dev_y_respect_to_c(m_val,c_val,k_val,t_vec[i],y0_val,y0_dot_val)
        
        #Evaluating f derivative respect to c  
        sk_num=Dev_y_respect_to_k(m_val,c_val,k_val,t_vec[i],y0_val,y0_dot_val)
        
        sens_vec[i,0]=sc_num
        
        sens_vec[i,1]=sk_num
       
    res=np.matmul(np.transpose(sens_vec),sens_vec) #Sensitivy matrix [2]
    
    return res

"""
Comentarios para explicar algunas de las funciones definidas en este archivo:
    
-Este archivo es auxiliar. Al importarlo desde un archivo principal permite
utilizar todas las funciones que hayan sido definidas en este archivo.

-Sensitivity_matrix2 (...) es la función que permite construir la matriz de 
sensitividad cuando los coeficientes a determinar son c y k (en ese orden). Es,
por lo tanto, una matriz de 2x2. Puede verificarse al ejecutar este archivo
si primero se hacen "activas" las lineas del código desde la línea 134 en adelante. 
Para ello se selecciona todas esas lineas y en la ventana Edit se selecciona
"Uncomment" para que las líneas se conviertan en instrucciones o código.


-En las diagonales de la matriz de sensitividad están las estimaciones que se
requieren para calcular las varianzas de c y de k. Así, la raíz cuadrada
de la entrada de la matriz de sensitividad ubicada en 0,0  es el requerido
para estimar un intervalo de confiabilidad para un valor de c estimado.

-De forma similar, la raíz cuadrada de la entrada de la matriz de sensitividad
ubicada en 1,1 es el requerido para estimar un intervalo de confiabilidad para 
un valor de k estimado.

-Al obtener estos valores se procede de forma similar al ejemplo Example2_2nn.py
para estimar los intervalos de confiabilidad para c y k, respectivamente.
"""


  


    
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


# sen_mat=Sensitivity_matrix2(m_val,c_val,k_val,t_vec,y0_val,y0_dot_val)




   


   
