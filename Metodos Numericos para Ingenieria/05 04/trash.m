clc
%fun = @(x) sin (x)
fun = @(x) (x.^2)
x0 = 3
x = fzero(fun,x0)