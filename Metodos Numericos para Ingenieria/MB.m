function c=MB(a,b,tol,mI)
  c=((a+b)/2)
  i=0
  while i<mI & abs(f(c))>tol;
    if f(a)*f(c)<0;
      b=c;
    else
      a=c;
    endif
    c=(a+b)/2;
    i=i+1;
  endwhile
  c
endfunction