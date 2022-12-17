function [w,amp] = Fourier(t,y)
% calcula la velocidad angular w y la amplitud amp de una onda dada por t e
% y mediante la Transformada Radipa de Fourier
N=length(t);
w0=2*pi/N; 
real(1:N+1)=0; 
im(1:N+1)=0; 
dt = t(2) - t(1);
for k=0:N-1 
    for n=0:N-1 
        k1=k+1; 
        n1=n+1; 
        angle=k*w0*n;
        real(k1)=real(k1)+y(n1)*cos(angle)/N; 
        im(k1)=im(k1)-y(n1)*sin(angle)/N;
    end 
    amp(k1)=sqrt(real(k1)^2+im(k1)^2); 
    w(k1)=(k*w0/dt)/(2*pi);
end