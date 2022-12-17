function [x,e] = SRS (A,b,X_0,Tol,MaxIter)
  n=length(b)
  x(n,1)
  w=(2/(1+sen(pi/(n-1))))
  k=1
  while  k<=MaxIter;
    k++;
    for i=1:n;
      x(i)=(w/A(i,i))*(b(i)-A(i,1:i-1)*x(1:i-1) - A(i,i+1:n)*x(i+1:n))+(1-w)*x(i)
    endfor
    e(k)=(norm(x-x_0,inf)/norm(x,inf))
    if e(k)>Tol;
      x_0=x
  endwhile 
endfunction