function [w,amp] = MFourier(t,y)
%Calcula las frecuencias y amplitudes en base a la propia de MatLab para
%Tranformada rapida de Fourier
Y = fft(y);
N = length(y);
P = abs(Y/N);
amp = P(1:end);
amp(2:end-1) = 2*amp(2:end-1);
Fs = N/t(end);
w = Fs*(0:N-1)/N;
w = w';