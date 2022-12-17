#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 16:20:08 2022

@author: nekko

Ejemplo 1 entregado por el profesor para trabajr en clase
"""
"""
IDMI310-22
Introducción a la cuantificación de incertidumbres en Ingeniería
Magister en Ingeniería Mecánica y Materiales
Facultad de ciencias de la ingeniería
Universidad Austral de Chile
________________________________________
Created on Mon Mar 28 11:33:01 2022


Example 1 [1]:
A 1 DOF system to illustrate a simple example of unicertainty propagation, an
example taken from [1] and use it according to GNU General Public License 
(see license a the UQ project file).

If you want to review Python basics you may check ref. [2] or another reference.


    
References:
[1] Bilionis,  Ilias . "Hands-On Activity 1.1: The Uncertainty Propagation Problem",
ME 539 - Introduction to Scientific Machine Learning at
https://github.com/PredictiveScienceLab/data-analytics-se#me-539----introduction-to-scientific-machine-learning    

[2] Bilionis,  Ilias . Preface — Introduction to Data Science for Mechanical Engineers (Lecture Book).”
https://purduemechanicalengineering.github.io/me-297-intro-to-data-science/index.html (accessed Sep. 28, 2022).


Assume that following information is known [1]:
    
k:trailer's equivalent stiffness coefficient----> Manufacturing uncertainty [159,999, 160,001] N/m

v:trailer speed---->Operating condition: [80, 150] km/hour

m: trailer's equivalent mass---->Loading condition:[100, 200] kg

road surface is modelled as a sinusoidal
y_0: rough road's amplitude---->Road condition: [0, 100] mm     

L: wavelength ---->Road condition: [1, 2] m


It would be assumed that each input variable is a random variable distributed uniformly.

One way to generate samples of a x random variable such as x~U([a,b)])(it means x is uniformly
distributed over the interval [a,b)] ), is based on a Z random variable such as Z~U([0,1)]) is:
                                                                       
     x= a+(b-a)Z 

"""

#Importing the required libraries  
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
#_______________________________________

# The number of samples we wish to take
num_samples = 1000000 #It would be the same result for higher values?

#Creating arrays to store the samples we will take

k_samples = np.ndarray((num_samples, ))  #k
m_samples= np.ndarray((num_samples, ))   #m

y0_samples= np.ndarray((num_samples, ))  #y_0
v_samples= np.ndarray((num_samples, ))   #v
lam_samples=np.ndarray((num_samples, ))  #L

omegas = np.ndarray((num_samples, ))

#Output (QoI: Quantity of interest)
Xs = np.ndarray((num_samples, ))    # To store the samples of Xs

#Sampling
for i in range(num_samples):
    k = 160000.0 + np.random.rand()                 # np.random.rand() samples a number uniformly between 0 and 1
    m = 100.0 + (200.0 - 100.0) * np.random.rand()  # Here we sample a random number in [100, 200]
    y0 = 0.0+ (100.0-0.0) * np.random.rand() * 1e-3              # Turning it to m
    v = (80.0 + (150.0 - 80.0) * np.random.rand()) * 1e3 / 3600.0 # Turning it to m/s
    lam = 1.0 + (2.0 - 1.0) * np.random.rand()   #[m]
    omega = 2.0 * np.pi * v / lam        #[rad/s]
    
    #Amplitude response for steady state condition [1]
    X = np.abs(k * y0 / (k - m * omega ** 2)) #[m]
    
    #Storing each sample:
    k_samples[i]=k
    m_samples[i]=m
    y0_samples[i]=y0
    v_samples[i]=v
    lam_samples[i]=lam

    omegas[i] = omega
    
    Xs[i] = X

# pdfs based on kernel density estimation

kernel = stats.gaussian_kde(Xs)
x_pts=np.linspace(Xs.min(),Xs.max(),100)
estimated_pdf=kernel.evaluate(x_pts)

##Plotting historgrams for each variable

# Plot the k samples
fig, ax = plt.subplots(dpi=150)
ax.hist(k_samples, density=True)
ax.set_xlabel('k [N/m]')
ax.set_ylabel('$p(k)$ (Probability density of k)')

# Plot the m samples
fig, ax = plt.subplots(dpi=150)
ax.hist(m_samples, density=True)
ax.set_xlabel('m [Kg]')
ax.set_ylabel('$p(m)$ (Probability density of m)')    

# Plot the y0 samples
fig, ax = plt.subplots(dpi=150)
ax.hist(y0_samples*1e3, density=True)
ax.set_xlabel('y_0 [mm]')
ax.set_ylabel('$p(y_0)$ (Probability density of $y_0$)')    

# Plot the v samples
fig, ax = plt.subplots(dpi=150)
ax.hist(v_samples*3600.0/1e3, bins=100, alpha=0.5, density=True)
ax.set_xlabel('v [Km/h]')
ax.set_ylabel('$p(v)$ (Probability density of v)')    

# Plot the lam samples
fig, ax = plt.subplots(dpi=150)
ax.hist(lam_samples, density=True)
ax.set_xlabel('L [m]')
ax.set_ylabel('$p(L)$ (Probability density of L)')    


# Plot the angular velocity
fig, ax = plt.subplots(dpi=150)
ax.hist(omegas, bins=100, alpha=0.5, density=True)
ax.set_xlabel('$\omega$ [rad/s]')
ax.set_ylabel('$p(\omega)$ (Probability density of $\omega$)')

# Plot the amplitude response
fig, ax = plt.subplots(dpi=150)
ax.hist(Xs, density=True)
ax.set_xlabel('$X$ (Amplitude)')
ax.set_ylabel('$p(X)$ (Probability density of amplitude)')    

#Plot the pdf
fig, ax = plt.subplots(dpi=150)
ax.plot(x_pts,estimated_pdf)
ax.set_xlabel('$X$ (Amplitude)')
ax.set_ylabel('$p(X)$ (Probability density of amplitude)')