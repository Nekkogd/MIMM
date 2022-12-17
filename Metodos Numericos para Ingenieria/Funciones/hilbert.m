function H=matrix(n)
  H=[];
  for i=1:n;
    for j=1:n;
      H(i,j)=1/(i+j-1);
    endfor
  endfor
  H
endfunction
