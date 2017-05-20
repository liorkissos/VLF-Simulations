
clc
clear
close all

%% OLD
fs = 200;
T = 1/fs;
t = 0:T:1;
 delta = 2;
 a = 0.1;
%% OLD


%%% NEW
% fs=250000;
% F_if=10e3;
% BW=10e3;
% T_chip=0.2;
% 
% t1=10*T_chip;
% f0=F_if-BW/2;
% f1=F_if+BW/2;
% %beta=BW/t_max; % assures that freq sweeps over [F_if-BW/2,F_if+BW/2]

% Ts=1/fs;
% t=0:Ts:(round(t1/Ts)-1)*Ts;
%%% NEW

%%% Signal anf matched filter
%x = chirp(t,5,1,20);
x=chirp(t,f0,t1,f1);
h = fliplr(x);

figure
set(gcf,'windowstyle','docked')

subplot(3,2,1)

plot(t,x);grid;
xlabel('Time (s)')
title('Chirp waveform, x(t)')


subplot(3,2,2);
plot(t,h);
grid;
xlabel('Time (s)')
title('Matched filter, h(t)')

N=6*fs;                     %  Signal is measured for up to 6 sec.
s = [zeros(1,delta*fs) a*x];
s = [s zeros(1,N-length(s))];
t1 = linspace(0,6,1200);
%noise = 2*std(a*x)*randn(size(s));
noise = 0*std(a*x)*randn(size(s));
y =s + noise;

subplot(3,2,3);
plot(t1,y); grid;% Received signal y(t) is completely buried in noise and so is not noticeable  
xlabel('Time (s)')
title('Received signal, y(t)')
maxlag = 5*fs;             % maximum lag is defined up to 5 seconds
[Rxy, tau] = xcorr(y,x,maxlag);
tau = tau/fs;

subplot(3,2,4);
plot(tau (maxlag +1 : end), Rxy(maxlag +1 : end));grid;
xlabel(' Lag (\tau) ')
title('Cross-correlation \itR_x_y\rm(\tau)')
title('Cross-correlation using xcorr');
ym = conv(y,h);  % ym(t) is calculated by performing convolution of y(t) and h(t).
ym = ym(1:length(y));

subplot(3,2,5);
plot(t1(1: maxlag),ym(1:maxlag)); grid;
xlabel('Time (s)')
title('Filtered signal \ity_m\rm(t)')

% f=linspace(-fs/2,fs/2,length(y));
% A=fftshift(fft(y));
% plot(f,db(abs(A)))

