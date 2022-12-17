function [f,i] = CW(r,Re,Tol,Imax)
## CW(r,Re)calcula el factor de friccion en funcion de r(rugosidad relatica) y Re(Numero de Reynolds)
  fa=64/Re
  y(x)=1/(-2*log10((r/3.7)+(2.51/(Re*sqrt(x))))
  #c = 0
  while c <= Imax;
    ra = fa - y(fa)/numderivative(y, fa);
    errorabsoluto = abs((ra-fa)/ra)*100;
    fa = ra;
    c = c + 1;
    if(errorabsoluto < Tol);
      break
    end
  end
  i = c
endfunction
