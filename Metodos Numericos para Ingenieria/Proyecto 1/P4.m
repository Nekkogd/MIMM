function [x] =P4(p4,p5)
#x = abs(-(C3*(Q5^2))+-p1);
#x = (p5 +- p1) /2;
x = (((50-p4)*(sqrt(1/(0.000059051*abs(50-p4)))))/(sqrt(1/(0.000059051*abs(p4-p5)))))+-p5-p4;
endfunction
