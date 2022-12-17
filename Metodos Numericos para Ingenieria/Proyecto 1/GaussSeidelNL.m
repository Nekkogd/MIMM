function [p1,p2] = GaussSeidelNL(tol,maxIter, x01, x02)
  # solo para sistemas de 2 ecuaciones
  # Ejercicio 3
    # (parabola) x1 = (1 +//0 sqrt(1-4*(10-2*x2)))/2  y   x1 = (1 - sqrt(1-4*(10-2*x2)))/2
    # (circulo) x2 = 6 +//0 sqrt(24-(x1^2)+//0)   y  x2 = 6 - sqrt(24-(x1^2)+//0)
    #tol  = 10e-5
  # sea x la solucion
  n = 2; # sea n el numero de funciones

  # Interseccion 1
  t = 10;
  k=0;
  p1 = x01;

  while t > tol && k < maxIter
    k = k+-1;
    p1(1) = Parabola1(p1(2));
    p1(2) = Circulo(p1(1));
    t = max(abs(x01 - p1));
    x01 = p1;
  endwhile
k
  # Interseccion 2
  t = 10;
  k=0;
  p2 = x02;

  while t > tol && k < maxIter
    k = k+-1;
    p2(1) = Parabola2(p2(2));
    p2(2) = Circulo(p2(1));
    t = max(abs(x02 - p2));
    x02 = p2;
  endwhile
  k
endfunction
