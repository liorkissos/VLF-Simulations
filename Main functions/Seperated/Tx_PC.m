
%function [ Bits_Error,Bits_compared_total,SNR_other ] = OFDM_Link_function(Coding_config,OFDM_config,IF_chain_config,Channel_config,Simulation_config)


%if nargin<1 %% User Input %%%%%%%%%%%%%%%%%%

%% Self called
%
close all
clear
clc

%test_signal_processing_flag=1;
test_signal_processing_flag=0;



%% USER Input
%%%%%%%%%%%%%%%%%%

%% Simulation parameters

Link_Type='HW'
%Link_Type='SW'

Configuration='Operational' % OFDM, real channel, non identical symbols, IF signal, Minn& Zeng time synchronization


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

%%% D/A ports
DAC_port=0;
DAC_FS=10; % Volt

%%% Control the D/A's Full scale voltage
Vpeak_out=[]; % default 10v
%Vpeak_out=0.01;


%%% compensate for the Coil's 1/jwL transfer function in the digital domain
Pre_emphasis_flag=1;
%Pre_emphasis_flag=0;

%% PHY Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~Ref_Cfg_flag
    
    %% OFDM Settings
    
        %%% Config # 1 : 512 suncarriers- for long delay spreads: CP
        %%% extremely long (N_FFT/4) and LS equalizer
%         F_chip=20e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
%         %   F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
%         N_FFT=512; % do not vary!
%         Npilots=6;
%         %Npilots=8;
%         Nguard_band_left=28*2; % do not vary!Nguard_band_right=Nguard_band_left-1; % do not vary!
%         Nguard_band_right=Nguard_band_left-1; % do not vary!
%         % N_CP=50; % ML eqaulzier
%         N_CP=128; % LS equalizer
%         Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
%         P_total=1; % total OFDM symbol (time domain) power
%         M=16; % QAM order
%         N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
%         N_preamble_synch=4; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
%         Enhancement_prmbl_CE=5.2;
%         Enhancement_prmbl_synch=4.1;
    
%     %%% Config # 4 : 802.16a
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
    M=64; % QAM order
   % M=16; % QAM order
    N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
    N_preamble_synch=8; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
    Enhancement_prmbl_CE=4.3;
    Enhancement_prmbl_synch=3.2;
    
    
    %%% Equalizer
    
    %Equalizer_type='LS'
    %Equalizer_type='ML'
    
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
    t_coding=10; % 64 QAM.number of assured correction (in messages per codeword terms)
    %t_coding=4; % 16QAM
    
    
else
    
    load('PHY_params_Ref_Cfg.mat')
    disp('reference system configuration')
    
end

%% Audio objects definitions


T_recorded=1; % required recording duration in seconds

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


if DAC_FS<9.2
    warning('DAC output level tool low')
end

if N_preamble_CE<2
    error('2 preamble_CE are needed for noise estimation')
end


if N_CP>N_FFT/4
    error('N_CP is too large. may create problem of time synchronization: may detect the peak created by the CP autocorrelation')
end


if N_FFT*N_preamble_synch<2048
    error('Preamble Synch (long preamble) too short. might create frequency and timing correction problems')
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
OFDM_config.Equalizer_type=[];
OFDM_config.N_preamble_synch=N_preamble_synch;
OFDM_config.Enhancement_prmbl_CE=Enhancement_prmbl_CE;
OFDM_config.Enhancement_prmbl_synch=Enhancement_prmbl_synch;

IF_chain_config.F_if=F_if;
IF_chain_config.Frec_req=Frec_req;
IF_chain_config.Fs_req=Fs_req;
IF_chain_config.Pre_emphasis=Pre_emphasis_flag;

Audio_config.Voice_flag=Voice_flag;
Audio_config.q_bits=q_bits;
Audio_config.Fs_Audio=Fs_Audio;

Simulation_config.Input_Type=Input_Type;
Simulation_config.Time_sync_flag=Time_sync_flag;
Simulation_config.BB_IF_flag=BB_IF_flag;
Simulation_config.Configuration=Configuration;
Simulation_config.Frequency_sync_flag=Frequency_sync_flag;
Simulation_config.N_symbols=N_symbols;
%Simulation_config.Design_flag=Design_flag;
Simulation_config.Ref_Cfg_flag=Ref_Cfg_flag;
Simulation_config.Identical_symbols=Identical_symbols;

%% Simulation
%%%%%%%%%%%%%%%%

%% 1) Generation & QAM Modulation

if Generation_flag
    Audio_config.Voice_flag=0;
    bit_stream=[];
    
    N_data_bits_sSymbol=N_data*log2(M); % data bit per super symbol (=OFDM symbol)
    LCM=lcm(N_data_bits_sSymbol,lcm(8,K_coding));
    bit_stream=randi([0 1],LCM,1); %
    
    
    [ Symbol_stream_Tx,data_PreCoded_Tx,data_Tx,OFDM_config,N_bits,t1_tx ] = Modulator(bit_stream, Coding_config,OFDM_config,Simulation_config,Audio_config);
else
    %load('Tx_stream_50k.mat')
  %  load('Tx_stream_10k_Nfft_256_Ncp_50_64QAM_tcode_10.mat')
   % load('Tx_stream_5k_Nfft_512_Ncp_128_16QAM_tcode_4.mat')
  load('Tx_stream_BW_10kHz_5k_Nfft_256_Ncp_50_64QAM_tcode_10.mat')
  %load('Tx_stream_BW_10kHz_5k_Nfft_256_Ncp_50_16QAM_tcode_4.mat')
 
    
end


%% Tranmission Loop

kk_max=2;
for kk=1:kk_max
    
    kk
    %%% 2) Transmission
    
    [Signal_Tx,Frec,OFDM_config,Testing_data,IF_chain_config] = Transmitter(Symbol_stream_Tx,Coding_config,OFDM_config,IF_chain_config,Simulation_config);
    %t2_tx=clock;

    %%% 3) D/A
    
    if   Frec>250e3
        error('Fsampling or Freconstruction is too high')
    end
   
    Fs=Fs_req; % the required Fsample is the actual Fsample.
    
    Signal_Tx_digital=Signal_Tx; % for analysis purposes later on
    
    %%% Termination frame creation: a random sequence lasting T_termination
    %%% seconds
    if kk==kk_max
        T_termination=0.1;
        N_termination=round(T_termination/(1/Frec));
        Signal_Tx= randi([0 1],N_termination,1); %;
        display('Termination Frame transmitted')
    end

    gain=DAC_HW( Signal_Tx,Frec,DAC_port,DAC_FS,Vpeak_out,F_if,F_chip);

    pause(3.5)
    
end

if test_signal_processing_flag
    
    Trec=1/Frec_req;
    figure
    set(gcf,'windowstyle','docked')
    t=0:Trec:Trec*(length(Signal_Tx)-1);
    plot(t/1e-3,Signal_Tx);xlabel('[msec]');grid on;grid minor
    title('Tx')
    xlabel('[msec]')
    
end
