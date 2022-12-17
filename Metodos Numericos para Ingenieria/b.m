function [b]=b(h,dx,n,Ta,T0,Tn)
  n=n/dx;
  b=(h*dx^2*Ta)*ones(n,1);
  b(1,1)= b(1,1)-T0;
  b(n,1)= b(n,1)-Tn;
endfunction