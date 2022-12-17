function [A] = mat3(h,dx,m)
  n=m/dx;
  r=((2+(h*dx^2))*ones(1,n));
  v=-1*ones(1,n-1);
  A=diag([r])+diag([v],-1)+diag([v],1);
endfunction