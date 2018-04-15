
clc
clear
close all

%% constants

c=3e8;

mu_0=4*pi*1e-7;
epsilon_0=8.8e-12;


Z_air=377
%Z_air=10

%% user defined

f=10e3;

%f=5e3:5e3:150e3;

%sigma=6e7; % copper

%%% conducting material properties
sigma=10e-3;

epsilon_r=7;
mu_r=1;


%% Intermediate calculations

w=2*pi*f;

epsilon=epsilon_0*epsilon_r;
mu=mu_0*mu_r;

lambda=c./f;



%% Calculations

Medium_conduction_metric= sigma./(w*epsilon); % a mteric of the conductivity of the material

if min(Medium_conduction_metric)>100
    disp('good conductor')
    delta_meters=sqrt(2./(w*mu*sigma))
else
    disp('poor conductor')
    delta_meters=(2/sigma)*sqrt(epsilon/mu)
end

Z_medium=sqrt(j*w*mu./(sigma+j*w*epsilon))

Transmission_dB=db(abs(2*Z_medium./(Z_medium+Z_air)),'voltage')

%% Display

if length(f)>1
    
    figure
    plot(f/1e3,delta)
    grid minor
    xlabel('Frequency [kHz]')
    ylabel('Skin depth [meters]')
    set(gcf,'windowstyle','docked')
    %title(['Skin effect Vs frequency: \sigma=',num2str(sigma/1e-3),'[mS/m]. \epsilon_r=',num2str(epsilon_r),''])
    %title({['Skin effect Vs frequency'],['  \sigma=',num2str(sigma/1e-3),'[mS/m]. \epsilon_r=',num2str(epsilon_r),''],['Conductivity metric at 10kHz: \sigma/\omega\epsilon=',num2str(min(Medium_conduction_metric)),'']})
    title({['Skin effect Vs frequency'],['  \sigma=',num2str(sigma/1e-3),'[mS/m]. \epsilon_r=',num2str(epsilon_r),''],['Conductivity metric at 10kHz: \sigma/\omega\epsilon=',num2str((Medium_conduction_metric(f==10e3))),'']})
    
else
    f_kHz=f/1e3
end


