
close all
clear
clc

dbstop if error

%dbstop if warning

clear global BBB

global BBB


%% USER Input
%%%%%%%%%%%%%%%%%%


%% Simulation parameters

BER_Thershold=1e-3; % if BW<20e3, lower to BER_Threshold=1e-3

Error_bits_threshold=400;
N_bits_Threshold=5e6;

%%% Simulation parameters

MIMO_depth=1;

Configuration='Operational' % OFDM, real channel, non identical symbols, IF signal, Minn& Zeng time synchronization
%Configuration='Calibration' % OFDM, 1-tap channel, identical symbols, BB signal, artificial time synchronization based on group delay summing along the chain and exact sampling times

N_symbols=10000; % number of symbols in the Frame

Ref_Cfg_flag=0;
%Ref_Cfg_flag=1;

%% Channel


%%% SNR
%SNR_method='EbNo'; % EbNo of the data symbols only!
SNR_method='EsNo'; % True Baseband signal SNR within signal's
%bandwidth: the SNR we would measure with a spectrum analyzer. i.e;
%before upsampling or IF upconversion (see Channel_SW function). EVM
%will eventually be different from EsNo value and
%will depend on length of signal.


%%% Multipath pattern

Delay_resp_vec=[1]; % 0=the time delay at which the first sample enters the channel
Amp_resp_vec=[1];
Phase_resp_vec=[0];


% Delay_resp_vec=[0 2]; %
% Amp_resp_vec=[1 0.5];
% Phase_resp_vec=[0 0];

% Delay_resp_vec=[0 4 ]; %
% Amp_resp_vec=[1 0.5 ];
% Phase_resp_vec=[0 0 ];

%% PHY Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~Ref_Cfg_flag
    
    %% OFDM Settings
    
    
    
    %%% Config # 4 : 512 suncarriers
    %         F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
    %         %   F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
    %         N_FFT=256*2 % do not vary!
    %         Npilots=6;
    %         %Npilots=8;
    %         Nguard_band_left=28*2; % do not vary!Nguard_band_right=Nguard_band_left-1; % do not vary!
    %         Nguard_band_right=Nguard_band_left-1; % do not vary!
    %         % N_CP=64;
    %         N_CP=50*1;
    %         Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
    %         P_total=1; % total OFDM symbol (time domain) power
    %         M=64; % QAM order
    %         N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
    %         N_preamble_synch=8; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
    %         Enhancement_prmbl_CE=5.2;
    %         Enhancement_prmbl_synch=4.2;
    
    
    
    %         %%% Config # 4 : 802.16a. 256 subcarriers
    F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
    %F_chip=2.94e3;
    N_FFT=256; % do not vary!
    Npilots=6;
    Nguard_band_left=28; % do not vary!Nguard_band_right=Nguard_band_left-1; % do not vary!
    Nguard_band_right=Nguard_band_left-1; % do not vary!
    % N_CP=64;
    N_CP=50;
    Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
    P_total=1; % total OFDM symbol (time domain) power
    M=64; % QAM order
    N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
    N_preamble_synch=8; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
    Enhancement_prmbl_CE=4.2; % coded
    Enhancement_prmbl_synch=3.3; % coded
%     Enhancement_prmbl_CE=4.0; % uncoded
%     Enhancement_prmbl_synch=3.0; % uncoded
     
    
    
    
    %%% Config # 3: 802.11a. 64 subcarriers
    %             F_chip=2.94e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)
    %             N_FFT=64;
    %             Npilots=2;
    %             Nguard_band_left=7; %  do not vary! lower number than 6 will harm the anti aliasing filter and thus the performance
    %             Nguard_band_right=Nguard_band_left-1; % do not vary!
    %             N_CP=12;
    %             Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
    %             P_total=1; % total OFDM symbol (time domain) power
    %             M=64; % QAM order
    %             N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
    %             N_preamble_synch=32; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
    %             Enhancement_prmbl_CE=7.15; %
    %             Enhancement_prmbl_synch=4.1; %
    
    
    
    
    %%% Equalizer
    
    %Equalizer_type='LS'
    Equalizer_type='ML'
    
    %% IF conversion & Interpolation (% BER performance is aligned with theory
    %%% up until F_if=16e3 (F_chip=1e4))
    %%%%%
    
    %%% configuration# 1
    %     F_if=14e3;
    %     Frec_req=250e3;
    %     Fs_req=100e3;
    
    
    %%% configuration # 2
        F_if=10e3;
        Frec_req=100e3;
        Fs_req=100e3;
    

    
    %% Coding
    Coding_flag=1;
    %Coding_flag=0
    
    Interleave_flag=1;
    %Interleave_flag=0;
    
    m_coding=M; % must be equal to the modulation depth, M. see comm.RSEncoder Help
    t_coding=10; % number of assured correction (in messages per codeword terms)
    
    % N_coding = 2^(log2(m_coding))-1;           % RS codeword length
    % K_coding = N_coding-2*t_coding;
    
else
    
    load('PHY_params_Ref_Cfg.mat')
    
end

%% Carrier Offset
%Carrier_offset=60*(F_if/1e6) % Tx to Rx Carrier (F_if) offset. F_if/1e6 means 1ppm
Carrier_offset=0 % Tx to Rx Carrier (F_if) offset. F_if/1e6 means 1ppm

%% Calibration mode config

if strcmp(Configuration,'Calibration')
    
    
    F_chip=5e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
    N_FFT=256; % do not vary!
    %Npilots=5;
    Npilots=6;
    %Nguard_band_left=24; % do not vary!Nguard_band_right=Nguard_band_left-1; % do not vary!
    Nguard_band_left=28; % do not vary!Nguard_band_right=Nguard_band_left-1; % do not vary!
    Nguard_band_right=Nguard_band_left-1; % do not vary!
    N_CP=50;
    Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
    P_total=1; % total OFDM symbol (time domain) power
    M=64; % QAM order
    N_preamble_CE=2;
    N_preamble_synch=8;
    Enhancement_prmbl_CE =3.5; % does not matter in cse of calibration, since no equalizer
    Enhancement_prmbl_synch =2.5;  % does not matter in cse of calibration, since no equalizer
    
    Equalizer_type='LS'; % does not matter in cse of calibration, since no equalizer
    
    %%% IF conversion & Interpolation (% BER performance is aligned with theory
    %%% up until F_if=16e3 (F_chip=1e4))
    %%%%%
    
    %%% Config reference: do not vary!
    %F_if=14e3;
    F_if=10e3;
    N_upsample=14; % do not vary in all Matlab simulation, because fits exactly the group delays! needs to allow reconstruction of the IF carried signal without aliasing
    N_upsample_ZOH=11; % do not vary in all Matlab simulation, because fits exactly the group delays!.D/A. an odd number, in order to have an integer total delay: ZOH+reconstruction filter
    Frec_req=100e3;
    Fs_req=100e3;
    
    
    Delay_resp_vec=[0]; % 0=the time delay at which the first sample enters the channel
    Amp_resp_vec=[1];
    Phase_resp_vec=[0];
    
    Carrier_offset=0; % Tx to Rx Carrier (F_if) offset. F_if/1e6 means 1ppm
    
    %Coding_flag=1 %
    Coding_flag=0 %
    
    Interleave_flag=1;
    %Interleave_flag=0;
    
    m_coding=M;
    t_coding=10; % 64QAM, rate 0.6
    %t_coding=3; % 16 QAM, ra6e 0.6
    N_coding = 2^(log2(m_coding))-1;           %
    K_coding = N_coding-2*t_coding;
    
else % oprational
    N_upsample=[]; %
    N_upsample_ZOH=11; %
    
    
end

%% Simulation Preliminaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% OFDM basic signal parameters: before Tx manipulations (Interpolation, Inverse Sinc etc.)

Length_Frame=N_FFT+N_CP; % length of the frame to be serialized->upsampled->reconstructed (D/A'ed)

N_data=N_FFT-(Nguard_band_left+Nguard_band_right+Npilots+1);

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

%%% Coding
N_coding = 2^(log2(m_coding))-1;           % RS codeword length
K_coding = N_coding-2*t_coding;

%% Simulation configuration

switch Configuration
    case 'Operational'
        
        Input_Type='OFDM';
        Time_sync_flag=1;
        Frequency_sync_flag=1;
        % Frequency_sync_flag=0;
        Frequency_sync_Coarse_flag=1;
        %Frequency_sync_Coarse_flag=0;
        BB_IF_flag='IF';
        Identical_symbols=0;
        Design_flag='Load';
        
    case 'Calibration'
        
        Input_Type='OFDM';
        Time_sync_flag=0;
        Frequency_sync_flag=0;
        Frequency_sync_Coarse_flag=0;
        BB_IF_flag='IF';
        % BB_IF_flag='BB';
        % Identical_symbols=1;
        Identical_symbols=0;
        Design_flag='Design';
        
        
        
    case 'Impulse Response'
        
        Input_Type='Impulse';
        Time_sync_flag=0;
        Frequency_sync_flag=0;
        Frequency_sync_Coarse_flag=0;
        
        BB_IF_flag='IF';
        Identical_symbols=0;
        Design_flag='Load';
        
    case 'Filters Design'
        
        Input_Type='OFDM';
        BB_IF_flag='IF';
        Frequency_sync_flag=1;
        Frequency_sync_Coarse_flag=1;
        Identical_symbols=0;
        Design_flag='Design';
        
        
        
end

%% Sanity checks

if(mod(Npilots,2)~=0)
    error('N_pilots is odd. change to even') % to enable uniform scattering of pilot subarriers
end


if(mod(N_data+1,Npilots-1)~=0)
    error('N_data_sub_carriers should be an Integer division of Npilots-1') % to enable uniform scattering of pilot subarriers
end
%

if(N_CP>N_FFT)
    error('N_CP should be smaller than N_FFT')
end

if mod(log2(N_FFT),1)~=0 %
    error('N_FFT should be a power of 2')
end


if N_CP>N_FFT/4
    error('N_CP is too large. may create problem of time synchronization: may detect the peak created by the CP autocorrelation')
end

if (N_FFT~=512 && N_CP/N_data>0.2587 || N_FFT==512 && N_CP/N_data>0.1270 ) && strcmp(Equalizer_type,'ML')
    error('N_CP is too large for ML equalizer. may create problem of Equalization; LS problem not sufficiently overdetermined for ML. Change to LS equalizer or reduce CP length ')
end

if N_preamble_synch<8
    error('Preamble Synch (long preamble) too short. might create frequency and timing correction problems')
end

if strcmp(Configuration,'Calibration') && F_chip~=5e3
    error('On calibration mode F_chip must remain 5kHz, otherwise the transmitter might change it without control. the transmiter does it in order to keep an integer number of delay and other things')
end

if strcmp(Configuration,'Calibration') && MIMO_depth~=1
    error('Calibration can only be done in a SISO link, by definition')
    
end

%% Preparations

%%% SNR Vectors in EsNo terms

switch M
    case 4
        
        if Coding_flag && numel(Delay_resp_vec)==1 % coded & flat channel
            SNR_vec=[0,4,6,8,10,12,14,16,18]; %
        else
            SNR_vec=[0,4,6,8,10,12,14,16,18]; % 4-QAM, no MP, Calibration
        end
        
    case 64
        
        if Coding_flag && numel(Delay_resp_vec)==1 % coded & flat channel  
            if MIMO_depth>1
                SNR_vec=[0,4,8,12,13,13.5,14,14.5]; % RS(63,39)
            else % MIMO_depth=1
                SNR_vec=[0,4,8,12,16,17,18,18.5,19,19.25,19.5,19.75,20]; % RS(63,39)
            end
        else % uncoded or non flat channel
            SNR_vec=[0,4,8,12,16,18,20,22,23,24]; %
        end
        
    case 16
        if Coding_flag && numel(Delay_resp_vec)==1 % coded & flat channel
            SNR_vec=[0,4,8,9,10,11,11.75,12.5,13]; %
        else
            SNR_vec=[0,4,8,10,12,14,16,17,18]; %
        end
        
end


if strcmp(SNR_method,'EbNo') % translate to EbNo terms
    SNR_vec=SNR_vec-10*log10(log2(M));
end

%% Preparation of structures

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

OFDM_config.PTS.PTS=0; % PTS is neutralized in such simulation

IF_chain_config.F_if=F_if;
IF_chain_config.N_upsample=N_upsample;
IF_chain_config.N_upsample_ZOH=N_upsample_ZOH;
IF_chain_config.Carrier_offset=Carrier_offset;
IF_chain_config.Frequency_sync_Coarse_flag=Frequency_sync_Coarse_flag;
IF_chain_config.Frec_req=Frec_req;
IF_chain_config.Fs_req=Fs_req;
IF_chain_config.Pre_emphasis=0; % needed only on HW mode with real coils

Channel_config.SNR_method=SNR_method;
%Channel_config.SNR=SNR;
Channel_config.Delay_resp_vec=Delay_resp_vec; % 0=the time delay at which the first sample enters the channel
Channel_config.Amp_resp_vec=Amp_resp_vec;
Channel_config.Phase_resp_vec=Phase_resp_vec;
Channel_config.Fs_channel=F_chip;

Simulation_config.Input_Type=Input_Type;
Simulation_config.Time_sync_flag=Time_sync_flag;
Simulation_config.BB_IF_flag=BB_IF_flag;
Simulation_config.Configuration=Configuration;
Simulation_config.Frequency_sync_flag=Frequency_sync_flag;
Simulation_config.N_symbols=N_symbols;
Simulation_config.Design_flag=Design_flag;
Simulation_config.Identical_symbols=Identical_symbols;
Simulation_config.Ref_Cfg_flag=Ref_Cfg_flag;
Simulation_config.Voice_only=[];
Simulation_config.MIMO_depth=MIMO_depth;

%%% Dispaly preparation: updates at every simulation loop
h_figure=figure(1);
set(h_figure,'WindowStyle','Docked')
subplot(1,2,1)
h_axes_BER_EbNo=gca;
h_line_BER_Ref_EbNo=semilogy([0:0.5:15],NaN(length([0:0.5:15]),1), 'g:o');
hold on
pause(1)
h_line_BER_EbNo=semilogy([0:0.5:15],NaN(length([0:0.5:15]),1), ':o');
hold on
grid on
grid minor
xlabel('EbNo data [dB]')
ylabel('BER')
xlim([0,18])
ylim(gca,[0.99e-5 1])

subplot(1,2,2)
set(h_figure,'WindowStyle','Docked')
h_axes_BER_EsNo=gca;
h_line_BER_Ref_EsNo=semilogy([0:0.5:15],NaN(length([0:0.5:15]),1), 'g:o');
hold on
pause(1)
h_line_BER_EsNo=semilogy([0:0.5:15],NaN(length([0:0.5:15]),1), ':o');
hold on
grid on
grid minor
xlabel('EsNo [dB]')
ylabel('BER')
xlim([0,18])
ylim(gca,[0.99e-5 1])

%% Simulation Loop
%%%%%%%%%%%%%%%%%%%%
kk=1;
BER=100; % dummy. to enable the first iteration of the coming while loop
while (kk<=length(SNR_vec) && BER>BER_Thershold)
    
    % EbNo_data=EbNo_data_vec(kk);
    SNR=SNR_vec(kk);
    
    Channel_config.SNR=SNR;
    
    % Aladding Modem Loop: runs the modem at a given SNR
    % until sufficient number of errors is accumulated or sufficient number
    % of bits is accumulated
    Bits_error_total=0;
    N_bits=0;
    
    while  (Bits_error_total<Error_bits_threshold && N_bits<N_bits_Threshold)
        
        [ Bits_Error,Bits_compared_total,SNR_other ] = OFDM_Link_function(Coding_config,OFDM_config,IF_chain_config,Channel_config,Simulation_config);
        
        Bits_error_total=Bits_error_total+Bits_Error;
        
        N_bits=N_bits+Bits_compared_total;
        
    end
    
    if strcmp(SNR_method,'EbNo')
        EbNo_data_vec(kk)=SNR;
        EsNo_vec(kk)=SNR_other;
    else
        EbNo_data_vec(kk)=SNR_other;
        EsNo_vec(kk)=SNR;
    end
    
    BER=Bits_error_total/N_bits;
    BER_vec(kk)=BER;
    BER_Ref(kk) = berawgn(EbNo_data_vec(kk),'qam',max(4,M));
    
    datetime('now')
    disp(['EbNo=',num2str(EbNo_data_vec(kk)),'[dB]. EsNo=',num2str(EsNo_vec(kk)),'[dB]. BER=',num2str(BER_vec(kk)),'. # of bits=', num2str(N_bits),'. # of bit errors=' ,num2str(Bits_error_total),'.'])% # of frames=',num2str(N_Frames),'. # of frame errors=',num2str(Error_frames_total),'']);
    
    %disp(['EbNo=',num2str(EbNo_data_vec(kk)),'[dB]. EsNo=',num2str(EsNo_vec(kk)),'[dB]. BER=',num2str(BER_vec(kk)),'. FER=',num2str(FER_vec(kk)),'. # of bits=', num2str(N_bits),'. # of bit errors=' ,num2str(Error_bits_total),'. # of frames=',num2str(N_Frames),'. # of frame errors=',num2str(Error_frames_total),'']);
    
    
    if numel(Delay_resp_vec)>1
        Amp_resp_vec=sort(Amp_resp_vec,'descend');
        Channel_title=['MP, ',num2str(numel(Delay_resp_vec)),'[taps].LOS/NLOS strongest=',num2str(round(db(Amp_resp_vec(1)/Amp_resp_vec(2)))),'[dB] . Delay spread=',num2str(max(Delay_resp_vec)-min(Delay_resp_vec)),...
            '[Tchip]'];
    else
        Channel_title='Flat';
    end
    
    
    if Coding_flag
        if Interleave_flag
            Coding_title=['RS(',num2str(N_coding),',',num2str(K_coding),'). Interleaved'];
        else
            Coding_title=['RS(',num2str(N_coding),',',num2str(K_coding),'). Non-Interleaved'];
        end
    else
        Coding_title='None';
    end
    
    if strcmp(Configuration,'Calibration')
        Equalizer_type_title='irrelevant (Calibration)';
    else
        Equalizer_type_title=Equalizer_type;
    end
    
    set(h_line_BER_Ref_EbNo,'XData',EbNo_data_vec(1:kk),'YData',BER_Ref(1:kk))
    pause(1)
    set(h_line_BER_EbNo,'XData',EbNo_data_vec(1:kk),'YData',BER_vec(1:kk))
    pause(1)
    %title (h_axes_BER_EbNo,{['BER Vs. EbNo data'],[' Mod=',num2str(M),'-QAM. FFT=',num2str(N_FFT),'. Tchip=',num2str((1/F_chip)/1e-3),' [msec]'],[' Channel=',Channel_title,'. Carrier Offset=',num2str(Carrier_offset/(F_if/1e6)),'[ppm]. Equalizer=',Equalizer_type_title,''],['Coding=',Coding_title,''],[' Fif=',num2str(F_if/1e3),' [kHz]. Fs=',num2str(Fs_req/1e3),' [kHz]. Frec=',num2str(Frec_req/1e3),' [kHz]']})
    title (h_axes_BER_EbNo,{['BER Vs. EbNo data. MIMO 1x',num2str(MIMO_depth),''],[' Mod=',num2str(M),'-QAM. FFT=',num2str(N_FFT),'. Tchip=',num2str((1/F_chip)/1e-3),' [msec]'],[' Channel=',Channel_title,''],['Carrier Offset=',num2str(Carrier_offset/(F_if/1e6)),'[ppm]. Equalizer=',Equalizer_type_title,''],['Coding=',Coding_title,''],[' Fif=',num2str(F_if/1e3),' [kHz]. BW=',num2str(F_chip/1e3),' [kHz]. Fs=',num2str(Fs_req/1e3),' [kHz]. Frec=',num2str(Frec_req/1e3),' [kHz]']})
    legend(h_axes_BER_EbNo,['Theoretical: AWGN, UnCoded,Channel-known '],['Actual'],'Location','southwest')
    xlim(h_axes_BER_EbNo,[min(EbNo_data_vec) min(EbNo_data_vec)+max(SNR_vec)-min(SNR_vec)+0.05])
    
    set(h_line_BER_Ref_EsNo,'XData',EsNo_vec(1:kk),'YData',BER_Ref(1:kk))
    pause(1)
    set(h_line_BER_EsNo,'XData',EsNo_vec(1:kk),'YData',BER_vec(1:kk))
    pause(1)
    title (h_axes_BER_EsNo,{['BER Vs. EsNo. MIMO 1x',num2str(MIMO_depth),''],[' Mod=',num2str(M),'-QAM. FFT=',num2str(N_FFT),'. Tchip=',num2str((1/F_chip)/1e-3),' [msec]'],[' Channel=',Channel_title,''],['Carrier Offset=',num2str(Carrier_offset/(F_if/1e6)),'[ppm]. Equalizer=',Equalizer_type_title,''],['Coding=',Coding_title,''],[' Fif=',num2str(F_if/1e3),' [kHz]. BW=',num2str(F_chip/1e3),' [kHz]. Fs=',num2str(Fs_req/1e3),' [kHz]. Frec=',num2str(Frec_req/1e3),' [kHz]']})
    %title (h_axes_BER_EsNo,{['BER Vs. EsNo'],[' Mod=',num2str(M),'-QAM. FFT=',num2str(N_FFT),'. Tchip=',num2str((1/F_chip)/1e-3),' [msec]'],[' Channel=',Channel_title,'. Carrier Offset=',num2str(Carrier_offset/(F_if/1e6)),'[ppm]. Equalizer=',Equalizer_type_title,''],['Coding=',Coding_title,''],[' Fif=',num2str(F_if/1e3),' [kHz]. Fs=',num2str(Fs_req/1e3),' [kHz]. Frec=',num2str(Frec_req/1e3),' [kHz]']})
    legend(h_axes_BER_EsNo,['Theoretical: AWGN, UnCoded,Channel-known '],['Actual'],'Location','southwest')
    xlim(h_axes_BER_EsNo,[min(EsNo_vec) min(EsNo_vec)+max(SNR_vec)-min(SNR_vec)+0.01])
    
    
    kk=kk+1;
end

