%script Ejercicio 1
%Resolvemos el problema de valor inicial mediante Runge-Kutta-Fehlberg
clc
clear all
close all
xint = [0,5];
y0 = -1;
[h,xi,yi] = RungeKuttaFehlberg(0,5,-1,10e-10,0.5,0.001,@fp2p1e1);
[xj,yj] = ode45('fp2p1e1',xint,y0);
figure
plot(xi,yi,'b-.',xj,yj,'m*')
title(['Obtenemos la siguente funcion con un h de:',num2str(h)])
legend({'Rutina Propia','Rutina MATLAB'},'Location','northwest')