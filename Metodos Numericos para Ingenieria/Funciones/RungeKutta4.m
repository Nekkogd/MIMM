function [ti,wi] = RungeKutta4(f,a,b,N,alpha)
%Metodo de Runge-Kutta de 4 pasos
%[ti,wi] = RungeKutta4(f,a,b,N,alpha)
%Donde a,b son extremos, N es el numero de intervalos y alpha es la 
% condicion inicial y f la función 
%Donde ti y wi son los vectores solución [ti == x, wi == y]
h = (b-a)/N
t = a;
w = alpha;
ti(1,1) = a;
wi(1,1) = alpha;
for i = 1:N
    K_1 = h*f(t,w);
    K_2 = h*f(t+(h/2), w+(K_1/2));
    K_3 = h*f(t+(h/2), w+(K_2/2));
    K_4 = h*f(t+h,w+K_3);

    w = w + ((K_1 + 2*K_2 + 2*K_3+ K_4)/6);
    t = a + i*h;
    ti(i+1,1) = t;
    wi(i+1,1) = w;
end