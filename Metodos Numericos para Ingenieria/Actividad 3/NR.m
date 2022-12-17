function [ra, i] = NR(func, ai, Tol, Imax)
    aprox= ai;
    while 1
        ra = aprox - func(aprox)/numderivative(func, aprox);
        errorabsoluto = abs((ra-aprox)/ra)*100;
        aprox = ra;
        i = i + 1;
        if(errorabsoluto < Tol | i  == Imax ) then
        break
        end
    end
endfunction
