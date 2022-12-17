function [x] = P2(p2,p5)
#x = abs(p1-(c1*(q1^2)));
#x = abs(((((p2-p5)/sqrt(c5*abs(p2-p5)))+-((p1-p2/sqrt(c1*abs(p1-p2)))))*sqrt(c2*abs(p2-p3)))+-p3);
x = ((((50-p2)*(sqrt(1/(0.000373263*abs(50-p2)))))+-(((p5-p2)*(sqrt(1/(0.000059051*abs(p5-p2))))))/(sqrt(1/(0.000373263*abs(p2-0)))))+-0-p2);
endfunction
