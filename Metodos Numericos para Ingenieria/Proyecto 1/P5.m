function [x] = P5(p2,p4,p5)
#x = abs(p2-(C5*(Q5^2)));
#x = (p2 +- p4)/2;
x = (((p4-p5)*(sqrt(1/(0.000059051*abs(p4-p5)))))/(sqrt(1/(0.000059051*abs(p5-p2)))))+-p2-p5;
endfunction
