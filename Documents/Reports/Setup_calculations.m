
clear
clc
close all

%% User input
f0=10e3;
r=1; % distance from loop

mu0=4*pi*1e-7;
Lcoil=485e-6;
Zcoil=2*pi*f0*Lcoil; % loop impedance
D=0.5; % loop diameter
Nwind_coil=20; % Loop # of windings

PAPR=13; % of signal
Gamp=15; % V/V. Maximum is 20 V/V
Vin_amp_peak=10; % D/A's output (input to Amplifier)

Sensitivity_sensor=0.05e9; % V/Tesla



%% Calculations

%%% Tx Coil
Vin_amp_rms=Vin_amp_peak*10^(-PAPR/20);
Vout_amp_rms=Vin_amp_rms*Gamp;
Iout_amp_rms=Vout_amp_rms/Zcoil;

%%% Magnetic field at distance r 
Brms_Tesla=mu0*pi*((D/2)^2)*Iout_amp_rms*Nwind_coil/(2*pi*r^3) ;
Brms_mGauss=(Brms_Tesla*1e4)/1e-3

%%% Rx Coil
Vout_sensor_rms=Brms_Tesla*Sensitivity_sensor
Pout_sensor_dBm=10*log10((Vout_sensor_rms^2/(2*50))/1e-3)
Pnoise_dBm=-174+10*log10(10e3);
SNR_dB=Pout_sensor_dBm-Pnoise_dBm;

% Pnoise_watt=10^9*(Pnoise_dBm/10)*1e3;

%% Sanity checks

if Gamp>20
    warning('Ampifier gain too high: maximum 20V/V')
end

if Vin_amp_peak*Gamp>155
    warning('output voltage from amplifier too high')
end

if Iout_amp_rms>7
    warning('Loop current is too high')
end

if Brms_mGauss>4
    warning('Magnetic field too high')
end


