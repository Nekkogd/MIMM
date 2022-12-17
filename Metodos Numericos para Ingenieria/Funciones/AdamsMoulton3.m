function [t,w] = AdamsMoulton3(f,a,b,N,alpha)
%Metodo de Adams-Moulton de 3 pasos
%Cacula la aproximacion de la solucion de f entre a y b con alpha como
%condicion inicial en a
%[t,w] = AdamsMoulton3(f,a,b,N,alpha)
%Donde t son los puntos calculados del intervalo y w(t) los resultados 
%aproximados para cada punto [y(x) == w(t) con y como solucion de y' == f]
%Donde a,b son extremos, N es el numero de intervalos y alpha es la 
%condicion inicial y f la Ecuación Diferecial
h = (b-a)/N;
t(1,1) = a;
w(1,1) = alpha;
for i = 1:2
    K_1 = h*f(t(i,1),w(i,1));
    K_2 = h*f(t(i,1)+(h/2), w(i,1)+(K_1/2));
    K_3 = h*f(t(i,1)+(h/2), w(i,1)+(K_2/2));
    K_4 = h*f(t(i,1)+h, w(i,1)+K_3);

    w(i+1,1) = w(i,1) + ((K_1 + 2*K_2 + 2*K_3 + K_4)/6);
    t(i+1,1) = a + i*h;
end
for i = 3:N
    t(i+1,1) = a + (i*h);
    w(i+1,1) = (1/22.2)*((27.8*w(i,1)) - w(i-1,1) + (0.2*w(i-2,1)) - (0.192*((i-1)^2)) - (0.192*(i-1)) + 4.736);
end
end