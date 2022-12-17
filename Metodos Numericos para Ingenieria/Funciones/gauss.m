function [A_n,b_n]=gauss(A,b)
  n=length(A);
  A(:,n+1)=b;
  for i=1:n-1;
    for j=i+1:n;
      m=(A(i,j)/A(i,i));
      A(j,:)=A(j,:)-m*A(i,:);
    endfor
  endfor
  b_n=A(:,n+1)
  A(:,n+1);
  A_n=A
endfunction
