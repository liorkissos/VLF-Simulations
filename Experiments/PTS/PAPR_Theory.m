
clc
clear
close all

debug=0;


%% User defined parameters

N1=128; % number of iterations
N2=256; % number of iterations


x_dB=1:0.1:12;

%%

x=db2pow(x_dB);

Pr_N1=1-(1-exp(-x)).^N1;
Pr_N2=1-(1-exp(-x)).^N2;
%Pr2=1-(1-exp(-x/2)).^N;


%% Display

figure
set(gcf,'windowstyle','docked')
%semilogy(x_dB,Pr1)
semilogy(x_dB,Pr_N1,x_dB,Pr_N2)
title('Pr(PAPR>PAPR0)')
grid minor
xlabel('PAPR0')
ylabel('Pr(PAPR>PAPR0)')
legend(['N=',num2str(N1),'[bins] '],['N=',num2str(N2),'[bins] '])
    


