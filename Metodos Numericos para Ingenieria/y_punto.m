function yp = y_punto(t,y) 
%YPUNTO y_punto = f(t,y) para problemas de valor inicial 
yp=10*exp(-((t-2).^2)/(2*0.075^2))-0.6*y;