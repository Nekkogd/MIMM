function [x] = Q2(p2,p3,C2)
x = abs((p3-p2)/sqrt(C2*abs(p3-p2)));
endfunction