
clear
clc
close all

%% User input

%%% Signal
f0=10e3;
BW=10e3;
PAPR=12; % of signal

r=6 ;% distance from loop

%%% Tx coil

%%% FESP-5135 (large)
Lcoil=485e-6;
D=0.5; % loop diameter
Nwind_coil=20; % Loop # of windings

%%% FESP-5133 1330 (Small)
% Lcoil=6.5e-3;
% D=0.125; % loop diameter
% Nwind_coil=225; % Loop # of windings


%%% Transmiting devices
Gamp=50; % V/V. Maximum is 20 V/V
Vin_amp_peak=10; % D/A's output (input to Amplifier)

%%% Rx Sensor
Sensitivity_sensor=0.05e9; % V/Tesla
MFSD_noise_sensor_Tesla=8e-6*1e-9; % Magnetic field spectral density [Tesla/sqrt(Hz)]

%%%% A/D
ENoB=12.5;
Fs=100e3; % Hz
V_FS_adc=10;

%%% Atmospheric noise
PSD_noise_atmspr_dBm=-174+125; % [dBm/Hz] relies on curve from Gibson, figure 6-2 average of curves B and C at 10kHz


%%% constants
mu0=4*pi*1e-7;

%% Calculations

%%%%%%%%%%%%%
%%% Signal

%%% Tx Coil
Zcoil=2*pi*f0*Lcoil; % loop impedance

Vin_amp_rms=Vin_amp_peak*10^(-PAPR/20);
Vout_amp_rms=Vin_amp_rms*Gamp;
Iout_amp_rms=Vout_amp_rms/Zcoil;

%%% Magnetic field at distance r 
Brms_Tesla=mu0*pi*((D/2)^2)*Iout_amp_rms*Nwind_coil/(2*pi*r^3) ;
Brms_mGauss=(Brms_Tesla*1e4)/1e-3;

%%% Rx Coil
Vout_sensor_rms=Brms_Tesla*Sensitivity_sensor;
Pout_sensor_dBm=10*log10((Vout_sensor_rms^2/(50))/1e-3) 

%%%%%%%%%%%%%
%%% Noise

%%% Sensor

VSD_noise_sensor=MFSD_noise_sensor_Tesla*Sensitivity_sensor; %V/sqrt(Hz)
PSD_noise_sensor=VSD_noise_sensor^2; % V^2/Hz
Vsqr_noise_sensor=PSD_noise_sensor*BW; % V^2
P_noise_sensor_dBm=10*log10(Vsqr_noise_sensor/(50)/1e-3);


%%% A/D
DR_adc=6*ENoB+1.7;
P_FS_adc_dBm=10*log10(V_FS_adc^2/(2*50)/1e-3);
P_noise_adc_dBm=(P_FS_adc_dBm-DR_adc)-10*log10(Fs/BW);  % P_FS_adc_dBm-DR_adc is the noise integrated across the whole nyquist zone. we need to scale it to the signa'l bandwidth

%%% Atmospheric noise
P_noise_atmspr_dBm=PSD_noise_atmspr_dBm+10*log10(BW)



 SNR_dB=Pout_sensor_dBm-max([P_noise_sensor_dBm,P_noise_adc_dBm,P_noise_atmspr_dBm])


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
    warning('Magnetic field at receiver too high')
end


