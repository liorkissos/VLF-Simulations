
clc
clear
close all

%% constants

c=3e8;

mu_0=4*pi*1e-7;
epsilon_0=8.8e-12;


%Z0=377;

Z_air=10

%% user defined

f=5e3;

%f=5e3:5e3:150e3;

%sigma=5.7e7; % copper

sigma=10e-3;

epsilon_r=7;
mu_r=1;

%% Simulation

omega=2*pi*f;

epsilon=epsilon_0*epsilon_r;
mu=mu_0*mu_r;

lambda=c./f;
    
Medium_conduction_metric= sigma./(omega*epsilon); % a mteric of the conductivity of the material


delta=sqrt(2./(omega*mu*sigma));

Z_medium=sqrt(j*omega*mu./(sigma+j*omega*epsilon));

Relection_dB=db(abs(2*Z_medium./(Z_medium+Z_air)),'voltage');

%% Display

figure
plot(f/1e3,delta)
grid minor
xlabel('Frequency [kHz]')
ylabel('Skin depth [meters]')
set(gcf,'windowstyle','docked')
%title(['Skin effect Vs frequency: \sigma=',num2str(sigma/1e-3),'[mS/m]. \epsilon_r=',num2str(epsilon_r),''])
title({['Skin effect Vs frequency'],['  \sigma=',num2str(sigma/1e-3),'[mS/m]. \epsilon_r=',num2str(epsilon_r),''],['Conductivity metric: \sigma/\epsilon\mu=',num2str(min(Medium_conduction_metric)),'']})


