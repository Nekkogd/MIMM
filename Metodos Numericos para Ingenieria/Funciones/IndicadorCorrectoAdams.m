function [t,w] = IndicadorCorrectoAdams(f,a,b,N,alpha)
%Donde a,b son extremos, N es el numero de intervalos y alpha es la 
% condicion inicial y f la función
clc
h = (b-a)/N
t(1,1) = a;
w(1,1) = alpha;
for i = 1:3
    K_1 = h*f(t(i,1),w(i,1));
    K_2 = h*f(t(i,1)+(h/2), w(i,1)+(K_1/2));
    K_3 = h*f(t(i,1)+(h/2), w(i,1)+(K_2/2));
    K_4 = h*f(t(i,1)+h, w(i,1)+K_3);

    w(i+1,1) = w(i,1) + ((K_1 + 2*K_2 + 2*K_3 + K_4)/6)
    t(i+1,1) = a + i*h
end
for i = 4:N
    t(i+1,1) = a + i*h;
    w(i+1,1) = w(i,1) + h*((55*f(t(i,1),w(i,1)) - 59*f(t(i-1,1),w(i-1,1)) + 37*f(t(i-2,1),w(i-2,1)) - 9*f(t(i-3,1),w(i-3,1)))/24);
    w(i+1,1) = w(i,1) + h*((9*f(t(i+1,1),w(i+1,1)) + 19*f(t(i,1),w(i,1)) - 5*f(t(i-1,1),w(i-1,1)) + f(t(i-3,1),w(i-3,1)))/24);
end
end