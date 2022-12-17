function [ra,n] = friccion(r,Re,Tol,Imax)
  r
  Re
  Tol
  Imax
  fa = 64/Re
  n = 0
  f(fa)=fcw(fa,r,Re)
  while n <= Imax;
    ra = fa - f(fa)/ numerative(f,fa);
    ea = abs((ra-fa)/ra)*100
    fa = ra
    n = n+1
    if er < Tol;
      break
    endif
  endwhile
endfunction
