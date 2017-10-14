
clc
clear
close all

%N=1e3;

T=1;

f0=1;

Fs=100*f0;

Ts=1/Fs;

N=T/Ts;

t=0:Ts:(N-1)*Ts;

% x1=cos(2*pi*f0*t);
% x2=cos(2*pi*2*f0*t);
% x3=cos(2*pi*3*f0*t);

x1=exp(j*2*pi*f0*t);
x2=exp(j*2*pi*2*f0*t);
x3=exp(j*2*pi*3*f0*t);



x=1*x1+j*x2-j*x3;


figure
set(gcf,'windowstyle','docked')
subplot(2,1,1)
plot(real(x1))
hold on
plot(real(x2))
hold on
plot(real(x3))
grid on

legend('g1','g2','g3')

%subplot(2,1,2)
title('The sum')
plot(real(x))
grid on