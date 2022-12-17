function [x] = Bi(f,a,b,tol,mI)
  x=(a+b)/2;
  k=0
  while k<mI & abs(f(c))>tol;
    if f(a)*f(x)<0;
      b=x
    else
      a=x
    endif
    x=(a+b)/2
    k=k+1
  endwhile
endfunction