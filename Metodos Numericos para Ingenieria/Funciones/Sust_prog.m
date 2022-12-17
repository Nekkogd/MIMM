function [y] = Sust_prog(L,b)
  n=length(b)
  y=zeros(n,1)
  y(1)=b(1,1)/L(1,1)
  for i=2:n;
    y(i)=(1/L(i,i))*(b(i,1)-L(i,1:i-1)*y(1:i-1))
  endfor
  y
endfunction