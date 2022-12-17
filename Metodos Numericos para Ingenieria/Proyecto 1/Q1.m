function [x] = Q1(p2,p1,C1)
x = (p2-p1)/sqrt(C1*abs(p2-p1));
endfunction