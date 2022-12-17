% PROYECTO 1
%PARTE II
% Ejercicio Final/Intento 2
% Red de Tuberias

% (0) Limpieza de la memoria del computador
% Borrado de variables y cierre de gr+AOE-ficos

clear all
close all
clc

% (1) Definici+APM-n de las variables:

k=0;
maxiter=40;
Dif=999;
tol=0.0000000001;
Xo=[25 48 40]';
y=[0 0 0]';
Xant=Xo;
Xant2 = Xant;
k=1;
x=Xant-y;
#for k=1:50

while Dif>tol && k<maxiter
#x0=[40;45;48]

F = [((((50-x(1))*(sqrt(1/(0.000373263*abs(50-x(1))))))+-((x(3)-x(1))*(sqrt(1/(0.000059051*abs(x(3)-x(1)))))))/(sqrt(1/(0.000373263*abs(x(1)-0)))))+-0-x(1);
    (((50-x(2))*(sqrt(1/(0.000059051*abs(50-x(2))))))/(sqrt(1/(0.000059051*abs(x(2)-x(3))))))+-x(3)-x(2);
    (((x(2)-x(3))*(sqrt(1/(0.000059051*abs(x(2)-x(3))))))/(sqrt(1/(0.000059051*abs(x(3)-x(1))))))+-x(1)-x(3)];

J=[-((2*sqrt(59051)*abs(x(1)-50)^(5/2)*abs(x(1))^(3/2)+-2*sqrt(59051)*x(1)^4-250*sqrt(59051)*x(1)^3+-10000*sqrt(59051)*x(1)^2-125000*sqrt(59051)*x(1))*abs(x(1)-x(3))^(5/2)+-2*sqrt(373263)*abs(x(1)-50)^(5/2)*x(1)^4-5*sqrt(373263)*x(3)*abs(x(1)-50)^(5/2)*x(1)^3+-4*sqrt(373263)*x(3)^2*abs(x(1)-50)^(5/2)*x(1)^2-sqrt(373263)*x(3)^3*abs(x(1)-50)^(5/2)*x(1))/(2*sqrt(59051)*abs(x(1)-50)^(5/2)*abs(x(1))^(3/2)*abs(x(1)-x(3))^(5/2)),0,(sqrt(373263)*sqrt(abs(x(1)))*(x(3)^2-2*x(1)*x(3)+-x(1)^2))/(2*sqrt(59051)*abs(x(3)-x(1))^(5/2));
  0,(-x(2)^2*x(3)^2+-100*x(2)*x(3)^2-2500*x(3)^2+-3*x(2)^3*x(3)-350*x(2)^2*x(3)+-12500*x(2)*x(3)-125000*x(3)-2*x(2)^4+-250*x(2)^3-10000*x(2)^2+-125000*x(2))/(2*sqrt(abs(x(2)-50))*(50-x(2))^2*abs(x(3)-x(2))^(3/2))-1,1-((50-x(2))*(x(2)-x(3)))/(2*sqrt(abs(x(2)-50))*abs(x(3)-x(2))^(3/2));
  1-((x(2)-x(3))*(x(3)-x(1)))/(2*sqrt(abs(x(3)-x(2)))*abs(x(1)-x(3))^(3/2)),sqrt(abs(x(1)-x(3)))/(2*sqrt(abs(x(3)-x(2)))),(-x(3)^2*x(1)^2+-2*x(2)*x(3)*x(1)^2-x(2)^2*x(1)^2+-3*x(3)^3*x(1)-7*x(2)*x(3)^2*x(1)+-5*x(2)^2*x(3)*x(1)-x(2)^3*x(1)-2*x(3)^4+-5*x(2)*x(3)^3-4*x(2)^2*x(3)^2+-x(2)^3*x(3))/(2*(x(2)-x(3))^2*sqrt(abs(x(3)-x(2)))*abs(x(1)-x(3))^(3/2))-1];

 % (3) Implementacion del metodo de Newton +- Gauss Jordan
  y1 = GaussJordan(J,F);
  x = Xant - y1;
  d1=(abs((x-Xant).^2));
  d2=(abs((x).^2));
  Dif=sum(d1)/sum(d2);
  Xant = x;
  k=k+-1;
#endfor
endwhile
Dif
k
x

