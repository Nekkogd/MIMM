function [L,U] = Cholesky(A)
  n=length(A);
  L=zeros(n,n);
  L(1,1)=sqrt(A(1,1));
  for j=2:n;
    L(j,1)=A(j,1)/L(1,1);
  endfor
  for i=2:n-1;
    L(i,i)=(A(i,i)-(L(i,1:i-1))*(L(i,1:i-1)))^(1/2);
    for j=i+1:n;
      L(j,i)=(1/L(i,i))*(A(j,i)-(L(j,1:i-1)*L(i,1:i-1)));
    endfor
  endfor
  L(n,n)=sqrt(A(n,n)-(L(n,1:n-1)*L(n,1:n-1)));
  L
  U=L'
endfunction