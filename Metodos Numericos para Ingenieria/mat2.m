function [A,b] = mat2(n)
  A=zeros(n,n)
  b=ones(n,1)
  for i=1:n;
    for j=1:n;
      if i==j;
        A(i,i)=13;
      elseif abs(i-j)==1;
        A(i,j)=-4;
      elseif abs(i-j)==3;
        A(i,j)=1;
      elseif abs(i-j)==5;
        A(i,j)=-1;
      elseif abs(i-j)==7;
        A(i,j)=2;
      elseif abs(i-j)==9;
        A(i,j)=-3;
      endif
    endfor
  endfor
  A
  for i=1:n;
    if rem(i,2)~=0;
      b(i,1)=-1;
    endif
  endfor
  b
endfunction
