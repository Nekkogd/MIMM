clear all
clc

x = linspace(1,10,10)
y = [1.3,3.5,4.2,5,7,8.8,10.1,12.5,13,15.6]
 
n = length(x)

s1 = 0
s2 = 0
s3 = 0
s4 = 0
  
for k =1:n
  s1 = s1 + y(k)*x(k)
  s2 = s2 + y(k)
  s3 = s3 + x(k)
  s4 = s4 + (x(k)^2)
endfor

c1 = ((n*s1)-(s2*s3))/((n*s4)-(s3^2))
c2 = ((s2*s4)-(s1*s3))/((n*s4)-(s3^2))

f = c1*x + c2

mccoef = polyfit(x,y,1)
f2 = mccoef(1)*x + mccoef(2)

#plot(x,y,'*r',x,f(:),'b')
plot(x,y,'*r',x,f,'b',x,f2,'g')