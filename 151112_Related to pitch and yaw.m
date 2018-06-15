clear all

R=5;
X=50000;
Theta=0:0.1:5;


Dy=R*tan(Theta*pi/180)+X-X*cos(Theta*pi/180);

Dx=R-R*cos(Theta*pi/180)+X*tan(Theta*pi/180);

plot(Theta, Dy-Dx);