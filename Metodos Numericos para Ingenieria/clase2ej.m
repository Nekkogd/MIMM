clear all
close all

a = -1;
b = 5;

f=@(x) (x.*cos(x)+(x.^3+1).*exp(-x));
int_teo = integral(f,a,b)