function [x] = GaussSeidel(A,b,x0,tol,maxIter)
clc
# sea A matriz de ecuaciones y x0 el vector inicial
# sea x la solucion

# Ejercicio 2
  #A = [4,1,-1,1;1,4,-1,-1;-1,-1,5,1;1,-1,1,3]
  #b = [-2;-1;0;1]
  #x0 = [0;0;0;0]
  #tol = 0.00001
  #maxIter = 100
  
#

n = length(x0);,
#tol = tol*ones(n,1);
t = 10
k=0;
x = x0;

s = 0;

  while t > tol && k < maxIter
    for i = 1:n
      for j = (i+1):n
        s = s + A(i,j)*x(j);
      endfor
      for j = 1:(i-1)
        s = s + A(i,j)*x(j);
      endfor
      x(i) = ((b(i) - s)/A(i,i));
      s = 0;
    endfor
    t = max(abs(x0 - x));
    k = k+1
    x0 = x
  endwhile
endfunction