% PROYECTO 1
% Ejercicio 4
% Metodo Newton Rhapson
% Ecuaciones No lineales

% (0) Limpieza de la memoria del computador
% Borrado de variables y cierre de gráficos

clear all
close all
clc


% (1) Definición de las variables:

k=0;
maxiter=100;
Dif=999;
tol=0.000001;
Xo=[-1 -2 1]';
y=[0 0 0]';
Xant=Xo;
x=Xant-y;
Xant2 = Xant
while Dif>tol && k<maxiter
  k = k+1;
  F=[((x(1)^3)+(x(1)^2)*x(2)-x(1)*x(3)+6)
  ((e^(x(1)))+(e^(x(2)))-(x(3)))
  ((x(2)^2)-2*x(1)*x(3)-4)];

  J=[(3*(x(1)^2)+2*x(1)*x(2)-x(3)) (x(1)^2) -x(1)
  (e^(x(1))) (e^(x(2))) -1
  -2*x(3) 2*x(2) -2*x(1)];

  y1 = GaussJordan(J,F);
  y2 = J\F;
  x = Xant - y1;
  x2 = Xant2 -y2;
  Dif = max(abs(x - Xant));
  Xant = x;
  Xant2 = x2;
endwhile
x
x2
k
