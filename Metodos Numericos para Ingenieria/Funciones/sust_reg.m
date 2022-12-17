function [x] = Sust_reg(U,y)
  n=length(y)
  x=zeros(n,1)
  x(n)=y(n,1)/U(n,n)
  for i=n-1:-1:1;
    x(i)=(1/U(i,i))*(y(i,1)-U(i,i+1:n)*x(i+1:n))
  endfor
  x
endfunction