clc
clear
close all

A=fftshift(fft([zeros(3,1);0.5;zeros(2,1);-0.75;zeros(3,1);1;zeros(4,1);0.25;zeros(2,1);-0.125],1e3));

plot(linspace(-1,1,length(A)),db(abs(A)))

xlabel('Radial Frequency [rad/sec/\pi]')
ylabel('Amplitude [dB]')
grid minor