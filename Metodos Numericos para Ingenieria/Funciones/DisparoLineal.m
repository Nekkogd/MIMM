function [y,yp] = DisparoLineal(a,b,alpha,beta,h,p,q,r)
%[y,yp] = DisparoLineal(a,b,alpha,beta,h,p,q,r)
% donde a y b son los extremos del interbalo, alpha y beta son las
% condiciones de frontera en a y b y N es el numero de subintervalos 
%aproxima la solucion de -y'' + p(x)y' + q(x)y +r(x) = 0

N = (b-a)/h;
u(1,1) = alpha;
u(2,1) = 0;
v(1,1) = 0;
v(2,1) = 1;

for i = 1:N
    x = a + (i-1)*h;
    k(1,1) =  h*u(2,i);
    k(1,2) = h*(p(x)*u(2,i) + q(x)*u(1,i) + r(x));
    k(2,1) = h*(u(2,i) + (1/2)*k(1,2));
    k(2,2) = h*(p(x+h/2)*(u(2,i) + (1/2)*k(1,2)) + q(x + h/2)*(u(1,i) + (1/2)*k(1,1)) + r(x + h/2));
    k(3,1) = h*(u(2,i) + (1/2)*k(2,2));
    k(3,2) = h*(p(x+h/2)*(u(2,i) + (1/2)*k(2,2)) + q(x + h/2)*(u(1,i) + (1/2)*k(2,1)) + r(x + h/2));
    k(4,1) = h*(u(2,i) + k(3,2));
    k(4,2) = h*(p(x+h)*(u(2,i) + (1/2)*k(3,2)) + q(x + h)*(u(1,i) + (1/2)*k(3,1)) + r(x + h));

    u(1,i+1) = u(1,i) + (1/6)*(k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1));
    u(2,i+1) = u(2,i) + (1/6)*(k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2));

    kp(1,1) =  h*v(2,i);
    kp(1,2) = h*(p(x)*v(2,i) + q(x)*v(1,i) + r(x));
    kp(2,1) = h*(v(2,i) + (1/2)*kp(1,2));
    kp(2,2) = h*(p(x+h/2)*(v(2,i) + (1/2)*kp(1,2)) + q(x + h/2)*(v(1,i) + (1/2)*kp(1,1)) + r(x + h/2));
    kp(3,1) = h*(v(2,i) + (1/2)*kp(2,2));
    kp(3,2) = h*(p(x+h/2)*(v(2,i) + (1/2)*kp(2,2)) + q(x + h/2)*(v(1,i) + (1/2)*kp(2,1)) + r(x + h/2));
    kp(4,1) = h*(v(2,i) + kp(3,2));
    kp(4,2) = h*(p(x+h)*(v(2,i) + (1/2)*kp(3,2)) + q(x + h)*(v(1,i) + (1/2)*kp(3,1)) + r(x + h));

    v(1,i+1) = v(1,i) + (1/6)*(kp(1,1) + 2*kp(2,1) + 2*kp(3,1) + kp(4,1));
    v(2,i+1) = v(2,i) + (1/6)*(kp(1,2) + 2*kp(2,2) + 2*kp(3,2) + kp(4,2));
end
w(1,1) = alpha;
w(2,1) = (beta - u(1,N+1))/v(1,N+1);

for i = 1:N
    w(1,i+1) = u(1,i+1) + w(2,1)*v(1,i+1);
    w(2,i+1) = u(2,i+1) + w(2,1)*v(2,i+1);
end
y = w(1,:);
yp = w(2,:);