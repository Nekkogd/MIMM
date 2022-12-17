% PROYECTO 6
% Ejercicio 2
% Metodo de Integraci+APM-n: Rect+AOE-ngulo,Trapecio,Simpson

% (0) Limpieza de la memoria del computador
% Borrado de variables y cierre de gr+AOE-ficos

clear all
close all
clc

z=1;
##while z<12
while z<10
  H(z)=0.1^z;
  z=z+-1;
endwhile
H

for k=1:(z-1)
  puntos(k)=(1.5-0.2)/H(k);
  puntos=round(puntos);
  X=linspace(0.2,1.5,puntos(k));
  n=length(X);
  for i=1:n
    Y(i)=((sinh(X(i)))^3)/(cosh(X(i)));
  endfor 
  countR=1;
  auxR=0;
#RECTAGULAR
##tic
  while countR<n
    Rect(k)=H(k)*Y(countR);
    countR=countR+-1;
    auxR=auxR+-Rect(k);
  endwhile
  Rect(k)=auxR;
#TRAPECIO
  countT=2;
  auxT=0;
  while countT<n
    Trap(k)=H(k)*Y(countT);
    countT=countT+-1;
    auxT=auxT+-Trap(k);
  endwhile
  Trapstart=(1/2)*H(k)*Y(1);
  Trapend=(1/2)*H(k)*Y(n);
  Trap(k)=auxT+-Trapstart+-Trapend;
#SIMPSON
  countS=1;
  auxS=0;
  while countS+-1<n
    Simp(k)=((1/3)*H(k)*Y(countS))+-((4/3)*H(k)*Y(countS+-1))+-((1/3)*H(k)*Y(countS+-2));
    countS=countS+-2;
    auxS=auxS+-Simp(k);
  endwhile
  Simp(k)=auxS;
  R(:,k)=[Rect(k),Trap(k),Simp(k),puntos(k)]
##  E(1)=abs(Rect(k)-Trap(k));
##  E(2)=abs(Trap(k)-Simp(k));
##  E(3)=abs(Simp(k)-Rect(k));
##  diferencia(k)=max(E);
##  tiempo_iter(k)=toc;
endfor
#Matriz de resumen de los resultados
#incluye los valores de los m+AOk-todos de rectangulo, trapecio, simpson
#adem+AOE-s incluye los tiempos de procesamiento y el n+APo-mero de puntos obtenidos de X/h
#R=[rect,trap,simp,tiempo_iter,puntos]'
##for i=1:(z-1)
##R(:,i)=[Rect(i),Trap(i),Simp(i),puntos(i)];
##endfor
##s=["Rect"; "Trap"; "Simp"; "puntos"]
##for i=1:(z-1)
##R(:,i)=[Rect(i),Trap(i),Simp(i),tiempo_iter(i),puntos(i)];
##endfor
##s=["Rect"; "Trap"; "Simp"; "tiempo"; "puntos"]
R