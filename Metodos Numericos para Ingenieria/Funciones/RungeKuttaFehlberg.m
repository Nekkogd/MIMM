function [h,x,y] = RungeKuttaFehlberg(a,b,alpha,TOL,hmax,hmin,f)
%Metodo de Runge-Kutta de 4 pasos
%[h,x,y] = RungeKuttaFehlberg(a,b,alpha,TOL,hmax,hmin,f)
%Donde a,b son extremos, alpha es la condicion inicial, TOL es la 
% tolerancia, hmax y hmin son los tamaños de paso maximo y minimo y f es la
%Ecuacion Diferencial
%Donde h es el paso y x e y son los vectores que aproximan la solución de f

t = a;
w = alpha;
h = hmax;
FLAG = 1;
p = 0;
while FLAG == 1
    p = p+1;
    K_1 = h*f(t,w);
    K_2 = h*f(t + (1/4)*h, w + (1/4)*K_1);
    K_3 = h*f(t + (3/8)*h, w + (3/32)*K_1 + (9/32)*K_2);
    K_4 = h*f(t + (12/13)*h, w + (1932/2197)*K_1 - (7200/2197)*K_2 + (7296/2197)*K_3);
    K_5 = h*f(t + h, w + (439/216)*K_1 - 8*K_2 + (3680/513)*K_3 - (845/4104)*K_4);
    K_6 = h*f(t + (1/2)*h, w - (8/27)*K_1 + 2*K_2 - (3544/2565)*K_3 + (1859/4104)*K_4 - (11/40)*K_5);

    R = (1/h) * abs((1/360)*K_1 - (128/4275)*K_3 - (2197/75240)*K_4 + (1/50)*K_5 + (2/55)*K_6);

    if R <= TOL
        t = t+h;
        w = w + (25/216)*K_1 + (1408/2565)*K_3 + (2197/4104)*K_4 - (1/5)*K_5;
        x(p) = t;
        y(p) = w;
    end
    delta = 0.84* ((TOL/R)^(1/4));
    if delta <= 0.1
        h = 0.1*h;
    elseif delta >= 4
        h = 4*h;
    else
        h = delta*h;
    end
    if h>hmax
        h = hmax;
    end
    if t >= b
        FLAG = 0;
    elseif t + h > b
        h = b - t;
    elseif h < hmin
        FLAG = 0;
        disp('h minima excedida, Procedimiento completado sin éxito')
    end
end