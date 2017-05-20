


%% Self called

close all
clear
clc

%test_signal_processing_flag=1;
test_signal_processing_flag=0;
%
% dbstop if error
% %dbstop if warning
%
clear persistent nn
clear global CCC
clear global BBB
%
global BBB

%% USER Input
%%%%%%%%%%%%%%%%%%

%% Simulation parameters

MIMO_depth=3;

Link_Type='HW';
%Link_Type='SW'

Configuration='Operational' ;% OFDM, real channel, non identical symbols, IF signal, Minn& Zeng time synchronization
%Configuration='Calibration' % OFDM, 1-tap channel, identical symbols, BB signal, artificial time synchronization based on group delay summing along the chain and exact sampling times
%Configuration='Impulse Response'; % signal containing impulses at each

N_symbols=10000; % number of QAM symbols in the Frame
% N_symbols=45590; % number of QAM symbols in the Frame


%Voice_flag=1;
Voice_flag=0;

%Design_flag=1; % set OFDM and IF parameters
%  Design_flag=0; % load saved OFDM and IF configuration

%Ref_Cfg_flag=1;
Ref_Cfg_flag=0;


%Generation_flag=1;
Generation_flag=0;

%% Hardware


% DAC_port=0;
% DAC_FS=9.5; % Volt
ADC_port=2;
ADC_FS=10;
ADC_Device='9215';
%ADC_Device='6212';

%% PHY Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~Ref_Cfg_flag
    
    %% OFDM Settings
    
        %%% Config # 1 : 512 suncarriers
%                 F_chip=20e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
%                 %   F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
%                 N_FFT=512; % do not vary!
%                 Npilots=6;
%                 %Npilots=8;
%                 Nguard_band_left=28*2; % do not vary!Nguard_band_right=Nguard_band_left-1; % do not vary!
%                 Nguard_band_right=Nguard_band_left-1; % do not vary!
%                 N_CP=128;
%                 %N_CP=50*1;
%                 Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
%                 P_total=1; % total OFDM symbol (time domain) power
%                 M=16; % QAM order
%                 N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
%                 N_preamble_synch=4; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
%                 Enhancement_prmbl_CE=5.2;
%                 Enhancement_prmbl_synch=4.1;
    
    %%% Config # 2 : 802.16a
    F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
    %   F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
    N_FFT=256; % do not vary!
    Npilots=6;
    %Npilots=8;
    Nguard_band_left=28; % do not vary!Nguard_band_right=Nguard_band_left-1; % do not vary!
    Nguard_band_right=Nguard_band_left-1; % do not vary!
    %N_CP=64;
    N_CP=50;
    Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
    P_total=1; % total OFDM symbol (time domain) power
   % M=64; % QAM order
    M=16; % QAM order
    N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
    N_preamble_synch=8; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
    Enhancement_prmbl_CE=4.3;
    Enhancement_prmbl_synch=3.2;
    
    
    %%% Equalizer
    
    %Equalizer_type='LS';
    Equalizer_type='ML';
    
    
    %% IF conversion & Interpolation
    
    %%% configuration# 1
    %     F_if=15e3;
    %     Frec_req=250e3;
    %     Fs_req=100e3;
    %
    %%% configuration # 2
    F_if=10e3;
    Frec_req=100e3;
    Fs_req=100e3;
    
    %         F_if=15e3;
    %         Frec_req=100e3;
    %         Fs_req=100e3;
    %
    %% Coding
    Coding_flag=1;
    %Coding_flag=0
    
    Interleave_flag=1;
    %Interleave_flag=0;
    
    m_coding=M; % must be equal to the modulation depth, M. see comm.RSEncoder Help
    %t_coding=10; % 64QAM. number of assured correction (in messages per codeword terms)
    t_coding=4; % 16QAM. number of assured correction (in messages per codeword terms)

    % burst_error=t_coding+30; % see Vishwantan, p.52
    
    
else
    
    load('PHY_params_Ref_Cfg.mat')
    disp('reference system configuration')
    
end


%% Carrier Offset

%Carrier_offset=60*(F_if/1e6) % Tx to Rx Carrier (F_if) offset. F_if/1e6 means 1ppm
Carrier_offset=0; % Tx to Rx Carrier (F_if) offset. F_if/1e6 means 1ppm

%% Audio objects definitions

Fs_Audio=44.1e3; % minimum rate to sample audio well is 8kHz

q_bits=16;% quantization level

%% Simulation Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OFDM basic signal parameters: before Tx manipulations (Interpolation, Inverse Sinc etc.)

Length_Frame=N_FFT+N_CP; % length of the frame to be serialized->upsampled->reconstructed (D/A'ed)

N_data=N_FFT-(Nguard_band_left+Nguard_band_right+Npilots+1); % DC subcarrier is a null too- that's why we reduce anotther data subcarrier

% T is the duration of the rectangular function of every g_n(t) basis
% function. Therefore the bandwidth of each subcarrier (a sinc function) is 1/T.
% The F_chip is the reconstruction frequency at which the OFDM modulator operatres after
% the S/P block that follows the IFFT&CP block. since the OFDM signal after
% (and also before) concatenation of the CP spans from -N_FFT/2T to +N_FFT/2T, the F_chip needs
% to be equal or bigger than 2*N_FFT/2T=N_FFT/T=> F_chip=N_FFT/T
% both F_chip and N_FFT are set by the user, thus T is derived from them
T=N_FFT/F_chip; % see explanation above: F_chip needs to N_FFT/T

T_chip=1/F_chip;% Sampling time of the OFDM transmitter/Reciver, before/after the IF Front-End

L=N_FFT+N_CP; % the number of chips per OFDM frame. chip: an output sample of the OFDM modulator to the Tx chain
T_frame=T_chip*L; % OFDM frame (OFDM symbol+CP duration) duration
T_QAM_symbol=T_frame/N_data; % QAM modulator symbol rate. needs to fit the chip rate. produces N_data_sub_carriers QAM symbols along T_frame[sec] (whereas L chips ae prduced along the same duration)

%%%  Bandwidth calculations: the signal in the frequency domain is composed
%%%  of concatenated N_FFT "Dinc" functions with a 1/T main lobe width and offset
%%%  by 1/T as well (Prof. Nazarathy lectures, p.25-32)

B_signal_Hz=(N_FFT/2)/T; % 1-sided signal's BW in Hz terms.
B_signal_net_Hz=(N_FFT/2-Nguard_band_right)/T; % Single sided (right side) actual signal BW In terms of Hz. Actual: without the guard band

%%%% Power constraints

Amp_pilots=db2pow(Amp_pilots_dB);

P_data=P_total/(Amp_pilots*Npilots+N_data);

P_preamble=P_total/floor(N_FFT);

%%% Time constants
T_payload_sym=T_chip*(N_FFT+N_CP);
T_preamble_CE=N_preamble_CE*T_payload_sym;
T_preamble_synch=T_chip*(N_preamble_synch*N_FFT+N_CP);

%%% Coding
N_coding = 2^(log2(m_coding))-1;           % RS codeword length
K_coding = N_coding-2*t_coding;           % RS message length

%% Simulation configuration


Input_Type='OFDM';
Time_sync_flag=1;
% Frequency_sync_flag=0;
Frequency_sync_flag=1;
Frequency_sync_Coarse_flag=1;
%Frequency_sync_Coarse_flag=0;
BB_IF_flag='IF';
Identical_symbols=0;
Design_flag='Load';

%% Sanity checks

%  if strcmp(Link_Type,'HW')

if strcmp(Configuration,'Calibration')
    error('Calibration is impossible on HW link: we cannot know the exact delay (needed in Calibration mode)')
end

if Fs_req~=100e3 || (Frec_req~=100e3 && Frec_req~=250e3)
    error('Fs or Frec have wrong values for HW simulation')
end

if F_chip==5e3 && Frec_req==250e3
    error('impossible to synchronize HW with BW=5kHz and Frec=250kHz')
end

%end

if(mod(Npilots,2)~=0)
    error('N_pilots is odd. change to even') % to enable uniform scattering of pilot subarriers
end


if(mod(N_data+1,Npilots-1)~=0) % including the DC
    error('N_data_sub_carriers should be an Integer division of Npilots-1') % to enable uniform scattering of pilot subarriers
end
%

if(N_CP>N_FFT)
    error('N_CP should be smaller than N_FFT')
end

if mod(log2(N_FFT),1)~=0 %
    error('N_FFT should be a power of 2')
end


if N_preamble_CE<2
    error('2 preamble_CE are needed for noise estimation')
end


if N_CP>N_FFT/4
    error('N_CP is too large. may create problem of time synchronization: may detect the peak created by the CP autocorrelation')
end

if (N_FFT~=512 && N_CP/N_data>0.2587 || N_FFT==512 && N_CP/N_data>0.1270 ) && strcmp(Equalizer_type,'ML')
    error('N_CP is too large for ML equalizer. may create problem of Equalization; LS problem not sufficiently overdetermined for ML. Change to LS equalizer or reduce CP length ')
end


if N_FFT*N_preamble_synch<2048
    error('Preamble Synch (long preamble) too short. might create frequency and timing correction problems')
end

if strcmp(Configuration,'Calibration') && F_chip~=5e3
    error('On calibration mode F_chip must remain 5kHz, otherwise the transmitter might change it without control. the transmiter does it in order to keep an integer number of delay and other things')
end
%% Preparations of structures

%%% Configuration Structures

Coding_config.Coding_flag=Coding_flag;
Coding_config.N_coding=N_coding;
Coding_config.K_coding=K_coding;
Coding_config.Interleave_flag=Interleave_flag;

OFDM_config.F_chip=F_chip;
OFDM_config.N_FFT=N_FFT;
OFDM_config.N_data=N_data;% to be sent to functions later
OFDM_config.Npilots=Npilots;
OFDM_config.Amp_pilots_dB=Amp_pilots_dB;
OFDM_config.Nguard_band_left=Nguard_band_left;
OFDM_config.Nguard_band_right=Nguard_band_right;
OFDM_config.N_CP=N_CP;
OFDM_config.M= M;
OFDM_config.P_total=P_total;
OFDM_config.T=T;
OFDM_config.B_signal_net_Hz=B_signal_net_Hz;
OFDM_config.P_data=P_data;
OFDM_config.N_preamble_CE=N_preamble_CE;
OFDM_config.Equalizer_type=Equalizer_type;
OFDM_config.N_preamble_synch=N_preamble_synch;
OFDM_config.Enhancement_prmbl_CE=Enhancement_prmbl_CE;
OFDM_config.Enhancement_prmbl_synch=Enhancement_prmbl_synch;


IF_chain_config.F_if=F_if;
%IF_chain_config.N_upsample=N_upsample;
%IF_chain_config.N_upsample_ZOH=N_upsample_ZOH;
IF_chain_config.Carrier_offset=Carrier_offset;
IF_chain_config.Frequency_sync_Coarse_flag=Frequency_sync_Coarse_flag;
IF_chain_config.Frec_req=Frec_req;
IF_chain_config.Fs_req=Fs_req;

% Channel_config.SNR_method=SNR_method;
% Channel_config.SNR=SNR;
% Channel_config.Delay_resp_vec=Delay_resp_vec; % 0=the time delay at which the first sample enters the channel
% Channel_config.Amp_resp_vec=Amp_resp_vec;
% Channel_config.Phase_resp_vec=Phase_resp_vec;

Audio_config.Voice_flag=Voice_flag;
Audio_config.q_bits=q_bits;
Audio_config.Fs_Audio=Fs_Audio;


Simulation_config.Input_Type=Input_Type;
Simulation_config.Time_sync_flag=Time_sync_flag;
Simulation_config.BB_IF_flag=BB_IF_flag;
Simulation_config.Configuration=Configuration;
Simulation_config.Frequency_sync_flag=Frequency_sync_flag;
Simulation_config.N_symbols=N_symbols;
Simulation_config.Design_flag=Design_flag;
Simulation_config.Ref_Cfg_flag=Ref_Cfg_flag;
Simulation_config.Identical_symbols=Identical_symbols;
Simulation_config.MIMO_depth=MIMO_depth;


%else  %called by an external function
%% Simulation
%%%%%%%%%%%%%%%%

if strcmp(Design_flag,'Load')
    load('Filters.mat')
end


%% Preperations
if ~Generation_flag
    %load('Tx_stream_10k.mat')
   % load('Tx_stream_10k_Nfft_256_Ncp_50_64QAM_tcode_10.mat')
   % load('Tx_stream_5k_Nfft_512_Ncp_128_16QAM_tcode_4.mat')
   %load('Tx_stream_BW_10kHz_5k_Nfft_256_Ncp_50_64QAM_tcode_10.mat')
   load('Tx_stream_BW_10kHz_5k_Nfft_256_Ncp_50_16QAM_tcode_4.mat')
  %  t1_rx=clock;
    Testing_data=[];
    Simulation_config.Voice_only=0;
else
    Testing_data=[];
    Audio_config.Voice_flag=1;
    Simulation_config.Voice_only=1;
    data_PreCoded_Tx=[];
    data_Tx=[];
end

T_termination=0.1;
EOF=0;
kk=0;
hError=comm.ErrorRate;
Bits_Error_total=0;
Bits_compared_total=0;

%% The reception l
while 1
    %% A/D
    timeout=10;
    [Signal_Rx1,aiSession,EOF]=ADC_HW( Frec_req,Fs_req,ADC_Device,ADC_port,ADC_FS,F_if,F_chip,[],'Separated',T_termination );
    
    kk=kk+1
    if EOF
        display('Communication terminated')
        pause(2)
        break
    else
        Signal_Rx= Signal_Rx1;
    end
    %% Reception
    Signal_Rx_MIMO=Signal_Rx;
    [ Symbol_stream_Rx ] = Receiver(Signal_Rx_MIMO,Fs_req,OFDM_config,IF_chain_config,Simulation_config,Testing_data );
    %[ Symbol_stream_Rx] = Receiver(Signal_Rx,Fs_req,OFDM_config,IF_chain_config,Simulation_config,Testing_data );

    %% QAM demodulation
    
    [ data_DeCoded_Rx,data_Rx,t2_rx ] = DeModulator(Symbol_stream_Rx,Coding_config,OFDM_config,Simulation_config,Audio_config );
    
    %% BER calculation
    if Coding_config.Coding_flag
        Error_strct=step(hError,data_PreCoded_Tx,data_DeCoded_Rx); % we ommit from the comparison the chopped symbols (due to group delay)
    else
        Error_strct=step(hError,data_Tx,data_Rx); % we ommit from the comparison the chopped symbols: due to group delay an fgthe fact that matlab's filter function chops the end of the vector, and we connot do anything about that
    end
    
    Bits_Error=Error_strct(2);
    Bits_compared=Error_strct(3);
    
    Bits_Error_total=Bits_Error_total+Bits_Error;
    Bits_compared_total=Bits_compared_total+Bits_compared;
end
%% Analysis

BER=Bits_Error_total/Bits_compared_total

Pavg=mean(Signal_Rx(:,1:MIMO_depth).^2,1)/50;
Pavg_dBm=db(Pavg/1e-3,'power')

hEVM=comm.EVM('Normalization','Average reference signal power');
%EVM=step(hEVM,Symbol_stream_Tx(1:length(Symbol_stream_Rx)),Symbol_stream_Rx);
EVM=step(hEVM,Symbol_stream_Tx,Symbol_stream_Rx);
EVM_dB=db(EVM/100)

length_test=round(length(Symbol_stream_Tx)/4);

EVM_start=step(hEVM,Symbol_stream_Tx(1:length_test),Symbol_stream_Rx(1:length_test));
EVM_start_dB=db(EVM_start/100);

EVM_end=step(hEVM,Symbol_stream_Tx(end-2*length_test:end-length_test),Symbol_stream_Rx(end-2*length_test:end-length_test));
EVM_end_dB=db(EVM_end/100);

f=linspace(-Fs_req/2,Fs_req/2,length(Signal_Rx));
SIGNAL_RX=fftshift(fft(Signal_Rx,[],1),1);

figure
set(gcf,'windowstyle','docked')
% plot(f_dca/1e3,db(abs(A)))
% hold on
plot(f/1e3,db(abs(SIGNAL_RX(:,1:MIMO_depth))))
grid on
grid minor
xlabel('frequncy [kHz]]')
title(['Rx after A/D. Fchip=',num2str(F_chip/1e3),'[kHz]. Fsample=',num2str(Fs_req/1e3),'[kHz]'])
legend({'A/D port=0','A/D port=1','A/D port=2','A/D port=3'})
%legend(cellstr(num2str('N','N=%-d')))
xlim([-Fs_req/2/1e3,Fs_req/2/1e3])


Ts=1/Fs_req;
t_Rx=0:Ts:Ts*(length(Signal_Rx)-1);

figure
set(gcf,'windowstyle','docked')
plot(t_Rx/1e-3,Signal_Rx(:,1:MIMO_depth))
xlabel('[msec]');grid on;grid minor
title(['Rx Signal at AFE interface.T preamble synch=',num2str(T_preamble_synch/1e-3),'[msec]. T preamble CE=',num2str(T_preamble_CE/1e-3),'[msec].'])
legend({'A/D port=0','A/D port=1','A/D port=2','A/D port=3'})

scatterplot(Symbol_stream_Rx);
grid on
grid minor
title (gca,{['EVM=',num2str(EVM_dB),'[dB]'],['MIMO 1x',num2str(MIMO_depth),''],[' Mod=',num2str(M),'-QAM. FFT=',num2str(N_FFT),'. Tchip=',num2str((1/F_chip)/1e-3),' [msec]'],['Equalizer=',Equalizer_type,''],[' Fif=',num2str(F_if/1e3),' [kHz]. BW=',num2str(F_chip/1e3),' [kHz]. Fs=',num2str(Fs_req/1e3),' [kHz]. Frec=',num2str(Frec_req/1e3),' [kHz]']})
set(gcf,'windowstyle','docked');

Psignal_average=mean(Signal_Rx(t_Rx/1e-3>50).^2);
Psignal_max=max(Signal_Rx(t_Rx/1e-3>50).^2);

PAPR_dB=db(Psignal_max/Psignal_average,'power')
