function [c,a,b] = FourierB(y)
%Calcula Transformada rapida de Fourier (Burden)
%[c,a,b] = FourierB(m,p,y) 
%donde c es la parte compleja; a y b son la parte real; m es ;p es ;y es

N = length(y);
m = N/2;
p = log2(m);
M = m;
q = p;
zeta = exp((pi*1i)/m);

for j = 1:(2*m)
    c(j) = y(j);
end
for j = 1:M
    xi(j+1)  = zeta^j;
    xi(j+M+1) = -xi(j+1);
end
K = 0;
xi(1) = 1;

for L = 1:(p+1)
    while K < (2*m-1)
        for j = 1:M
            K = 0;
            for g = p+1:1
                K = K + k(g)*(2^g);
            end
            K_1 = 0;
            K_1 = K/(2^q);
            K_2 = 0;
            for r = q:p
                K_2 = K_2 + k(r)*(2^(p-r));
            end
            eta = c(K+1+(M*xi*K_2));
            c(K+1+M) = c(K+1) - eta;
            c(K+1) = c(K+1) + eta;
            K = K +1;
        end
        K = K + M;
    end
    K = 0;
    M = M/2;
    q = q-1;
end
while K <(2*m-1)
    K = 0;
    j = 0;
    for h = p:0
        K = K + k(h)*(2^h);
        j = j + k(h-p)*(2^h);
    end
    if j > K 
        c(j) = c(k);
    end
end
a(1) = c(1)/m;
a(m+1) = real(((exp(-1i*pi*m))*c(m+1))/m);
for j = 1:(m-1)
    a(j+1) = real(((exp(-1i*pi*j))*c(j+1))/m);
    b(j+1) = imag(((exp(-1i*pi*j))*c(j+1))/m);
end