function F = myfun(x)
  F = [3*x(1) - cos(x(2)*x(3)) - 1/2; x(1)^2 - 81*((x(2)+0.1)^2) + sin(x(3)) + 1.06 ; exp(-x(1)*x(2))+ 20*x(3)+(10*pi-3)/3];