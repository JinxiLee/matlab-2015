clear all

A1=-0.0001107;

A2=0.09029;

A3=0.2313;

SR=1000;

tmax=125;

time=0:1/SR:tmax;

voltage=A1*time.^2+A2*time+A3;

plot(time,voltage);