%Ejercicio 2
clc
clear all
close all
a = 0;
b = 2;
%h = 0.2;
%N = (b-a)/h;
N = 10;
y0 = 0.5;

[xb,yb] = AdamsBashforth4(@fejemplo1,a,b,N,y0)
[xm,ym] = AdamsMoulton3(@fejemplo1,a,b,N,y0)
%[xr,yr] = RungeKutta4(@fejemplo1,a,b,N,y0);

figure
%plot(xb,yb,'b-.',xm,ym,'m*',xr,yr,'g-o')
plot(xb,yb,'b-.',xm,ym,'m*')
%title(['Obtenemos la siguente funcion con un h de:',num2str(h)])
legend({'Adams-Bashford de 4 pasos','Adams-Moulton de 3 pasos'},'Location','northwest')