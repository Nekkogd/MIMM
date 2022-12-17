function [l,V] = Potencia (A,v0,tol,maxIter)
  clc
  #Datos ejemplo
    #  A = [11 6 -2;-2 18 1;-12 24 13]
    #  v0 = [1 1 1]'

  # Datos Proyecto
    # A = [10,9,7,5;9,10,8,6;7,8,10,7;5,6,7,5]
    # v0 = [1 2 1 5]'
    # tol = 1e-5
    # maxIter = 100

  n = length(A); # determina la dimension del sistema
  V(:,1) = v0; #carga el vector inicial
  error = 10;
  c = 0;

  #Calculo del primer autovector y autovalor
  while error > tol && c < maxIter
    c = c+-1; # contador iteraciones
    v = A*V(:,1); # autovector no unitario
    g = norm(v); # autovalor
    v = v/g; # autovector unitario
    error = max(abs(V(:,1) - v));
    V(:,1) = v; # asigna el autovector a un amtriz de autovectores
  endwhile
  #g = round(g); # solo para sistemas con autovalores enteros
  l(1) = g; # asigna el autovalor a un amtriz de autovectores

  I = diag(ones(n,1));
  r = (A - ((l(1))*I));

  # Calculo de los demas autovalores y autovectores
  for h = 2:n;
      r = r*(A - ((l(h-1))*I));
      V(:,h) = r*v0;
      c = 0;
      error = 10;
      while error > tol && c < maxIter
        c = c+-1;
        v = A*V(:,h);
        g = norm(v);
        v = v/g;
        error = max(abs(V(:,h) - v));
        V(:,h) = v;
      endwhile
      #g = round(g); # solo para sistemas con autovalores enteros
      l(h) = g;
  endfor
  #l # vector de autovalores
  #V # matriz de autovectores
endfunction
