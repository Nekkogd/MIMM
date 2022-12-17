function [x] = P52(p2,p3,p5,Q1,C2)
x = (p2-((((p2-p3)/sqrt(C2*abs(p2-p3)))-Q1)*(sqrt(C2*abs(p2-p5)))));
endfunction