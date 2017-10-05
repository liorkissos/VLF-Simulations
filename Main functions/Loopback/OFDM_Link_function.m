
function [ Bits_Error,Bits_compared_total,SNR_other ] = OFDM_Link_function(Coding_config,OFDM_config,IF_chain_config,Channel_config,Simulation_config)


if nargin<1 %% User Input %%%%%%%%%%%%%%%%%%
    
    %% Self called
    
    close all
    clear
    clc
    
    %test_signal_processing_flag=1;
    test_signal_processing_flag=0;
    
    dbstop if error
    %dbstop if warning
    
    clear persistent nn
    clear global BBB
    
    global BBB
    
    %% USER Input
    %%%%%%%%%%%%%%%%%%
    
    %% Simulation parameters
    
    MIMO_depth=1;
    
    
    %Link_Type='HW'
    Link_Type='SW'
    
    Configuration='Operational' % OFDM, real channel, non identical symbols, IF signal, Minn& Zeng time synchronization
    %Configuration='Calibration' % OFDM, 1-tap channel, identical symbols, BB signal, artificial time synchronization based on group delay summing along the chain and exact sampling times
    %Configuration='Impulse Response'; % signal containing impulses at each
    
    N_symbols=30000; % number of QAM symbols in the Frame
    % N_symbols=45590; % number of QAM symbols in the Frame
    
    %%% PTS
    PTS=1 % PTS enabling flag
%     M_PTS=16; %number of divisions of the N_FFT long block 
%     L_PTS=4; % upsamling rate of the PAPR teting during the algorithm execution
%     W_PTS=8; % number of phase factors. e.g; for W_PTS=4, it is 1,j,-1,-j
    
    M_PTS=8; %number of divisions of the N_FFT long block 
    L_PTS=4; % upsamling rate of the PAPR teting during the algorithm execution
    W_PTS=4; % number of phase factors. e.g; for W_PTS=4, it is 1,j,-1,-j
    
    PTS_algorithm= 'Reduced_Complexity';
    %PTS_algorithm= 'Iterative_Flipping';
    
    scrambling= 'contiguous';
    %scrambling='interleaved';
    
    %Voice_flag=1;
    Voice_flag=0;
    
    %Design_flag=1; % set OFDM and IF parameters
    %  Design_flag=0; % load saved OFDM and IF configuration
    
    %Ref_Cfg_flag=1;
    Ref_Cfg_flag=0;
    
    %  Load_channel_flag=1;
    Load_channel_flag=0;
    
    %% Hardware
    
    
    DAC_port=0;
    DAC_FS=10; % Volt
    ADC_port=0;
    ADC_FS=10;
    ADC_Device='9215';
    %ADC_Device='6212';
    
    %% PHY Configuration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if ~Ref_Cfg_flag
        
        %% OFDM Settings
        
%         %%% Config # 1 : 512 suncarriers- for long delay spreads: CP
%         %%% extremely long (N_FFT/4) and LS equalizer
%         F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
%         %   F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
%         N_FFT=512; % do not vary!
%         Npilots=6;
%         %Npilots=8;
%         Nguard_band_left=28*2; % do not vary!Nguard_band_right=Nguard_band_left-1; % do not vary!
%         Nguard_band_right=Nguard_band_left-1; % do not vary!
%         N_CP=50; % ML eqaulzier
%         %N_CP=128; % LS equalizer
%         Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
%         P_total=1; % total OFDM symbol (time domain) power
%         M=64; % QAM order
%         N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
%         N_preamble_synch=4; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
%         
%         
%         if ~PTS
%             Enhancement_prmbl_CE=5.2;
%             Enhancement_prmbl_synch=4.1;
%         else
%             Enhancement_prmbl_CE=2.7; %
%             Enhancement_prmbl_synch=2; %
%             
%         end

        
        
        
        
%                 %%% Config # 2 : 802.16a. 256 subcarriers
%                 F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)N_FFT=128;
%                 %F_chip=2.94e3; % can be interchanged with 20e3 or any other bandwidth as long as N_FFT remains 256
%                 N_FFT=256; % do not vary!
%                 Npilots=6;
%                 Nguard_band_left=28; % do not vary!Nguard_band_right=Nguard_band_left-1; % do not vary!
%                 Nguard_band_right=Nguard_band_left-1; % do not vary!
%                 %N_CP=64;
%                 N_CP=50;
%                 Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
%                 P_total=1; % total OFDM symbol (time domain) power
%                 M=64; % QAM order
%                 N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
%                 N_preamble_synch=8; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
%                 
%                 if ~PTS
%                     Enhancement_prmbl_CE=4.3; %
%                     Enhancement_prmbl_synch=3.2; %
%                 else
%                     Enhancement_prmbl_CE=2; %
%                     Enhancement_prmbl_synch=0.5; %
%                     
%                 end


        
        
        
                %%% Config # 3: 802.11a. 64 subcarriers. When enabling
                %%% PTS, use it and not the others since greater N_FFT demand greater splitting (M_PTS) which becomes
                %%% too complex to realize
                % F_chip=2.94e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)
                F_chip=10e+03;% The signal's sampling frequency at the output of the cP insertion block. (Original value is: 1.3913e+08)
                N_FFT=64;
                Npilots=2;
                Nguard_band_left=7; %  do not vary! lower number than 6 will harm the anti aliasing filter and thus the performance
                Nguard_band_right=Nguard_band_left-1; % do not vary!
                %N_CP=12;
                N_CP=12;
                Amp_pilots_dB=0; % pilot subcarrier power vs average data subcarrier power
                P_total=1; % total OFDM symbol (time domain) power
                M=64; % QAM order
                N_preamble_CE=2; % at least 2 are needed for SNR calculation in receiver
                N_preamble_synch=32; % the length of the time domain long preamble: (N_preamble_synch*N_FFT+N_CP)*T_chip. do not go below 8! needed at low SNR's
                
                if ~PTS
                    Enhancement_prmbl_CE=7.15; %
                    Enhancement_prmbl_synch=4.1; %
                else
                    Enhancement_prmbl_CE=3.5; %
                    Enhancement_prmbl_synch=0; %
                end
        

        
        
        %%% Equalizer
        
        %Equalizer_type='LS'
        Equalizer_type='ML'
        
        
        %% IF conversion & Interpolation
        
        %%% configuration# 1
        %         F_if=15e3;
        %         Frec_req=250e3;
        %         Fs_req=100e3;
        %
        %%% configuration # 2
        F_if=10e3;
        Frec_req=100e3;
        Fs_req=100e3;
        
        
        %%% configuration # 2
        %         F_if=10e3;
        %         Frec_req=44.1e3;
        %         Fs_req=44.1e3;
        %
        %% Coding
        Coding_flag=1
        %Coding_flag=0
        
        Interleave_flag=1;
        %Interleave_flag=0;
        
        m_coding=M; % must be equal to the modulation depth, M. see comm.RSEncoder Help
        t_coding=10; %  64QAM. number of assured correction (in messages per codeword terms)
        %t_coding=4; %  16QAM. number of assured correction (in messages per codeword terms)
        
        
        
    else
        
        load('PHY_params_Ref_Cfg.mat')
        disp('reference system configuration')
        
    end
    
    %% Channel
    %%% Software
    
    %%% SNR
    % Important: This is the Signal to Noise ratio calculated of the whole nyquist band!
    %i.e: from -Fs/2 to +Fs/2 (-Fnyq to +Fnyq). It is also called the SNRin, the SNR of the
    %received r(t) signal ( r(t)=s(t)+v(t)), BEFORE it is being filtered by the
    %receiver's filter. The actual SNR, after the filtering will be higher by
    % 10log10(Fnyq/BW of signal single sided)
    % settings of the AWGN block: 1) measure with the spectrum analyzer
    % (channel measurements) the signal's power in dBW (!) 2) set the AWGN
    % block "Input signal power( reference to 1 Ohm)"to db2pow(signal power
    % [dBW]) 3) set the SNR to
    % EbNo-10*log10(Fnyq_Air_Hz/(B_signal_Hz/2))+10*log10(log2(M)) , where EbNo
    % is set by the user
    
    %%% Important: what truly matters is whether (on calibration mode): 1)
    %%% on EbNo mode: EbNo~EVM-log10(log2(M)) (e.g; M=4&EbNo=30, EVM=-33 on every signal length)
    %%% 2) on EsNo_data mode: EsNo~EVM(e.g; EsNo_data=30, EVM=-30 on every constellation and every signal length)
    %%% 3) BER curve is well aligned with theory
    
    %SNR_method='EbNo'; % SNR per bit on data  bits only! same as Esno_data apart from the log10(log2(M)) factor
    SNR_method='EsNo_data'; % Calibration parameter. on Calibration mode, EVM~EsNo_data. identical to EbNo apart from being lower by log10(log2(M)) regarding EbNo for the same configuration. e.g (calibration mode); EsNo_data=30-> EVM~-30. EbNo=30-> EVM~-33
    %SNR_method='EsNo'; % True Baseband signal SNR within signal's
    %bandwidth: the SNR we would measure with a spectrum analyzer. i.e;
    %before upsampling or IF upconversion (see Channel_SW function). EVM
    %will eventually be different from EsNo value and
    %will depend on length of signal. why? because overhead symbols, the
    %preambles, count more when total stream length
    %is short. Normally, SNR per data symbol
    %(EsNo_data) is better than SNR per general symbol,
    %because the overhead symbols make us "waste"
    %signal energy without increasing data rate.
    %therefore , for reliable measurement, net dataastream
    %should be the longest possible in comparison with
    %overhead symbols length
    
    
    SNR=40;
    
    %%%% Multipath pattern:
    %%%%% Delay_resp_vec(k)=0: the time delay at which the 1st sample enters
    %%%%% the channel. a channel's tap that does not cause any delay to the
    %%%%% signal.
    %%%%% Delays are given Tchip quantities
    
    if ~Load_channel_flag % channel defined by user. taps are defined in T_chip terms
        
        Fs_channel=F_chip;
        
        Delay_resp_vec=[0 ]; %
        Amp_resp_vec=[1 ];
        Phase_resp_vec=[0 ];
        
        %         Delay_resp_vec=[0 10 15]; %
        %         Amp_resp_vec=[0.2 1 0.5 ];
        %         Phase_resp_vec=[0 0 0  ];
        %         %
        %                 Delay_resp_vec=[0 4 ]; %
        %                 Amp_resp_vec=[1 0.5 ];
        %                 Phase_resp_vec=[0 0 ];
        
        
    else % load channel model: Alon's room model
        
        Fs_channel=44.1e3;
        
        load('room_impulse_response')
        ind_begin=620; % cutting the impulse response: 1) CP is not long enough to be able to compensate for such a channel 2) simulation time gets longer as channel response gets longer
        ind_end=ind_begin+500; % derived channel's length needs to be shorter than CP; (ind_end-ind_begin)<(N_CP-offset_deliberate)*T_chip
        Amp_resp_vec=room_impulse_response(ind_begin:ind_end); % needs to be shorter than CP
        
        if test_signal_processing_flag
            Ts=1/Fs_channel;
            t=0:Ts:(length(room_impulse_response)-1)*Ts;
            t_derived1=t(ind_begin:ind_end)*F_chip;
            t_derived2=t(ind_begin:ind_end)/1e-3;
            
            figure
            subplot(3,1,1)
            set(gcf,'windowstyle','docked')
            plot(t/1e-3,room_impulse_response);grid on;grid minor
            xlabel('time [msec]')
            title('Channel impulse response- Full')
            subplot(3,1,2)
            set(gcf,'windowstyle','docked')
            plot(t_derived2,Amp_resp_vec);grid on;grid minor
            xlabel('time [msec]')
            title('Channel impulse response- Simulated')
            subplot(3,1,3)
            set(gcf,'windowstyle','docked')
            plot(t_derived1,Amp_resp_vec);grid on;grid minor
            xlabel('time [Tchip]')
            title(['Channel impulse response- Simulated. Check that shorter than CP!. N CP=',num2str(N_CP),'[Tchips]'])
            
            
            
            500*(1/Fs_channel)/(1/F_chip)
        end
        
        
        Delay_resp_vec=0:1:length(Amp_resp_vec)-1;
        Phase_resp_vec=zeros(1,length(Amp_resp_vec));
    end
    
    
    %%%% Carrier Offset
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
        Enhancement_prmbl_CE =4.3;
        Enhancement_prmbl_synch =3.4;
        
        Equalizer_type='LS'; % does not matter in cse of calibration, since no equalizer
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%% IF conversion & Interpolation (% BER performance is aligned with theory
        %%% up until F_if=16e3 (F_chip=1e4))
        %%%%%
        
        %%% Config reference: do not vary!
        F_if=14e3;
        N_upsample=14; % do not vary in all Matlab simulation, because fits exactly the group delays! needs to allow reconstruction of the IF carried signal without aliasing
        N_upsample_ZOH=11; % do not vary in all Matlab simulation, because fits exactly the group delays!.D/A. an odd number, in order to have an integer total delay: ZOH+reconstruction filter
        Frec_req=100e3;
        Fs_req=100e3;
        
        
        Delay_resp_vec=[0]; % 0=the time delay at which the first sample enters the channel
        Amp_resp_vec=[1];
        Phase_resp_vec=[0];
        
        Coding_flag=1; %
        m_coding=M; % must be equal to the modulation depth, M. see comm.RSEncoder Help
        
        t_coding=10; % 64QAM, rate 0.6
        % t_coding=3; % 16 QAM, ra6e 0.6
        % t_coding=3; % 16 QAM, ra6e 0.6
        N_coding = 2^(log2(m_coding))-1;           %
        K_coding = N_coding-2*t_coding;
        
        %Interleave_flag=1;
        Interleave_flag=0;
        
        
    else % oprational
        N_upsample=[]; %
        N_upsample_ZOH=11; %
        
        
    end
    
    %% Audio definitions
    
    %      %% OLD -start
    %
    %     Fs_Audio=8e3*1; % minimum rate to sample audio well is 8kHz
    %
    %     q_bits=8;% quantization level
    %
    %      %% OLD -end
    
    %%% NEW -start
    
    Fs_Audio=44.1e3; % minimum rate to sample audio well is 8kHz
    
    q_bits=16;% quantization level
    
    %%% NEW -end
    
    %% Simulation Preliminaries
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% OFDM basic signal parameters: before Tx manipulations (Interpolation, Inverse Sinc etc.)
    
    Length_Frame=N_FFT+N_CP; % length of the frame to be serialized->upsampled->reconstructed (D/A'ed)
    
    % if PTS
    %    N_data=N_FFT-(Nguard_band_left+Nguard_band_right+Npilots+M_PTS+1); % DC subcarrier is a null too- that's why we reduce anotther data subcarrier
    % else
    N_data=N_FFT-(Nguard_band_left+Nguard_band_right+Npilots+1); % DC subcarrier is a null too- that's why we reduce anotther data subcarrier
    % end
    
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
    
    
    %%% Rate
    bits_per_sSymbol=N_data*log2(M)*K_coding/N_coding
    Rate=bits_per_sSymbol/T_frame
    
    %% Simulation configuration
    
    switch Configuration
        case 'Operational'
            
            Input_Type='OFDM';
            Time_sync_flag=1;
            % Frequency_sync_flag=0;
            Frequency_sync_flag=1;
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
    
    if strcmp(Link_Type,'HW')
        
        if strcmp(Configuration,'Calibration')
            error('Calibration is impossible on HW link: we cannot know the exact delay (needed in Calibration mode)')
        end
        
        if Fs_req~=100e3 || (Frec_req~=100e3 && Frec_req~=250e3)
            error('Fs or Frec have wrong values for HW simulation')
        end
        
        if F_chip==5e3 && Frec_req==250e3
            error('impossible to synchronize HW with BW=5kHz and Frec=250kHz')
        end
        
    end
    
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
    
    if strcmp(Link_Type,'HW') && strcmp(Configuration,'Calibration')
        error('Calibration is impossible on HW link: we cannot know the exact delay (needed in Calibration mode)')
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
    
    if (N_FFT~=512 && N_CP/N_data>0.2587 || N_FFT==512 && N_CP/N_data>0.1270 ) && strcmp(Equalizer_type,'ML')
        error('N_CP is too large for ML equalizer. may create problem of Equalization; LS problem not sufficiently overdetermined for ML. Change to LS equalizer or reduce CP length ')
    end
    
    
    if N_FFT*N_preamble_synch<2048
        error('Preamble Synch (long preamble) too short. might create frequency and timing correction problems')
    end
    
    if strcmp(Configuration,'Calibration') && F_chip~=5e3
        error('On calibration mode F_chip must remain 5kHz, otherwise the transmitter might change it without control. the transmiter does it in order to keep an integer number of delay and other things')
    end
    
    if ~( M==64 && t_coding==10) && ~( M==16 && t_coding==4)
        error('wrong combination of modultion order and coding')
    end
    
    if N_FFT==64 && strcmp(Link_Type,'HW')
        error('N_FFT==64 on HW not working well for some reason. tried to debug but was unable to find resaon')
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
    
    OFDM_config.PTS.PTS=PTS;
    if PTS
        OFDM_config.PTS.M_PTS=M_PTS;
        OFDM_config.PTS.W_PTS=W_PTS;
        OFDM_config.PTS.L_PTS=L_PTS;
        OFDM_config.PTS.PTS_algorithm=PTS_algorithm;
        OFDM_config.PTS.scrambling=scrambling;
    end
    
    IF_chain_config.F_if=F_if;
    IF_chain_config.N_upsample=N_upsample;
    IF_chain_config.N_upsample_ZOH=N_upsample_ZOH;
    IF_chain_config.Carrier_offset=Carrier_offset;
    IF_chain_config.Frequency_sync_Coarse_flag=Frequency_sync_Coarse_flag;
    IF_chain_config.Frec_req=Frec_req;
    IF_chain_config.Fs_req=Fs_req;
    IF_chain_config.Pre_emphasis=0; % needed only on HW mode with real coils
    
    
    Channel_config.SNR_method=SNR_method;
    Channel_config.SNR=SNR;
    Channel_config.Delay_resp_vec=Delay_resp_vec; % 0=the time delay at which the first sample enters the channel
    Channel_config.Amp_resp_vec=Amp_resp_vec;
    Channel_config.Phase_resp_vec=Phase_resp_vec;
    Channel_config.Fs_channel=Fs_channel;
    
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
    Simulation_config.Voice_only=[];
    Simulation_config.MIMO_depth=MIMO_depth;
    
else  %called by an external function
    %% called by an external function
    
    Link_Type='SW';
    Voice_flag=0;
    test_signal_processing_flag=0;
    Audio_config.Voice_flag=0;
    
    F_chip=OFDM_config.F_chip;
    
    F_if=IF_chain_config.F_if;
    
    Fs_req=IF_chain_config.Fs_req;
    
    Design_flag=Simulation_config.Design_flag;
    
    MIMO_depth=Simulation_config.MIMO_depth;
    
    
end % if nargin <1

%% Simulation
%%%%%%%%%%%%%%%%

if strcmp(Design_flag,'Load')
    load('Filters.mat')
end

%% Generation & QAM Modulation

% if PTS
%     OFDM_config.N_data=N_data-M_PTS; %modulator is required to generate less than N_data in PTS case to leave space for the PTS coefficients
% end

bit_stream=[];

[Symbol_stream_Tx,data_PreCoded_Tx,data_Tx,OFDM_config,N_bits,t1_tx ] = Modulator(bit_stream,Coding_config,OFDM_config,Simulation_config,Audio_config);

% if PTS
%     OFDM_config.N_data=N_data; % we restore the original N_data
% end
%% Transmission


[Signal_Tx,Frec,OFDM_config,Testing_data,IF_chain_config] = Transmitter(Symbol_stream_Tx,Coding_config,OFDM_config,IF_chain_config,Simulation_config);
t2_tx=clock;


if strcmp(Link_Type,'HW')
    if Voice_flag
        T_AFE_iteration_estimated_sec=length(Signal_Tx)*(1/Frec); % based on emprical results; it takes the DAC always 0.2secs more than theory to empty its buffer
        Tspare=2;
        T_ADC_iteration_sec=T_AFE_iteration_estimated_sec+Tspare; % the ADC acquisition duration needs to be longer than DAC's, but shorter than the Audio device's
        T_Comm_iteration_estimated_sec=T_ADC_iteration_sec+0.45; % empirical results show that AFE HW need 0.45 sec in addition to assigned comm time period to terminate their action completely and to enable recording from Mic again
    else
        T_AFE_iteration_estimated_sec=length(Signal_Tx)*(1/Frec); % based on emprical results; it takes the DAC always 0.2secs more then theory to empty its buffer
        Tspare=0.3;
        T_ADC_iteration_sec=T_AFE_iteration_estimated_sec+Tspare; % the ADC acquisition duration needs to be longer than DAC's, but shorter than the Audio device's
    end
end

%% AFE & Wireless Channel propagation


switch Link_Type
    case 'SW'
        
        %%%% SISO- start
        %         Signal_Tx_digital=Signal_Tx; % for analysis purposes later on
        %         [Signal_Tx,Fs_channel,OFDM_config,Testing_data]=DAC_SW(Signal_Tx,Frec,IF_chain_config,OFDM_config,Simulation_config,Testing_data);
        %         [ Signal_Rx,SNR_other,Testing_data ] = Channel_SW( Signal_Tx,Fs_channel,Channel_config,Simulation_config,OFDM_config,Testing_data );
        %         [ Signal_Rx ,Fs,Testing_data] = ADC_SW( Signal_Rx,Fs_channel,IF_chain_config,Simulation_config,Testing_data );
        %         Signal_Rx_digital=Signal_Rx; % for analysis purposes later on
        %%%% SISO- end
        
        %%%% MIMO- start
        Signal_Tx_digital=Signal_Tx; % for analysis purposes later on
        [Signal_Tx,Fs_channel,OFDM_config,Testing_data]=DAC_SW(Signal_Tx,Frec,IF_chain_config,OFDM_config,Simulation_config,Testing_data);
        for mm=1:MIMO_depth
            [ Signal_Rx,SNR_other,Testing_data ] = Channel_SW( Signal_Tx,Fs_channel,Channel_config,Simulation_config,OFDM_config,Testing_data );
            [ Signal_Rx ,Fs,Testing_data] = ADC_SW( Signal_Rx,Fs_channel,IF_chain_config,Simulation_config,Testing_data );
            Signal_Rx_MIMO(:,mm)=Signal_Rx;
            Signal_Rx_digital=Signal_Rx; % for analysis purposes later on
            
        end
        %%%% MIMO- end
        
        
    case 'HW'
        
        %%% Timing calculations
        if Voice_flag
            T_AFE_iteration_estimated_sec=length(Signal_Tx)*(1/Frec); % based on emprical results; it takes the DAC always 0.2secs more than theory to empty its buffer
            Tspare=2;
            T_ADC_iteration_sec=T_AFE_iteration_estimated_sec+Tspare; % sent to HW.the ADC acquisition duration needs to be longer than DAC's, but shorter than the Audio device's
            T_Comm_iteration_estimated_sec=T_ADC_iteration_sec+0.45; % empirical results show that AFE HW need 0.45 sec in addition to assigned comm time period to terminate their action completely and to enable recording from Mic again
        else
            T_AFE_iteration_estimated_sec=length(Signal_Tx)*(1/Frec); % based on emprical results; it takes the DAC always 0.2secs more then theory to empty its buffer
            Tspare=0.3;
            T_ADC_iteration_sec=T_AFE_iteration_estimated_sec+Tspare; % sent to HW. the ADC acquisition duration needs to be longer than DAC's, but shorter than the Audio device's
        end
        
        
        Fs=Fs_req; % the required Fsample is the actual Fsample.
        
        if (strcmp(ADC_Device,'9215')&& Fs_req>100e3 ) || (strcmp(ADC_Device,'6212')&& Fs_req>250e3 ) || Frec>250e3
            error('Fsampling or Freconstruction is too high')
        end
        
        [ Aliasing_Margin] = Aliasing_check( Frec,Fs,F_if,1.2*(2*B_signal_Hz) ); % distance betwen closest replicas
        Signal_Tx_digital=Signal_Tx; % for analysis purposes later on
        
        T_ADC_iteration_sec=[]; % empty if you wan to synchronize with chirps
        [ Signal_Rx,Signal_Tx,gain ] = Channel_HW_w_synch( Signal_Tx,Frec,DAC_port,DAC_FS,Fs,ADC_Device,ADC_port,ADC_FS,F_if,F_chip,T_ADC_iteration_sec );
        Signal_Rx_MIMO=Signal_Rx;
        
        Signal_Rx_digital=Signal_Rx/gain; % for analysis purposes later on. Rx signal is subject to amplification applied on the Tx signal right before the transmission in order to exploit the full dynamic range of the D/A
end


Ts=1/Fs;

if test_signal_processing_flag
    figure
    subplot(2,1,1)
    t=0:Ts:Ts*(length(Signal_Tx)-1);
    plot(t/1e-3,Signal_Tx);xlabel('[msec]');grid on;grid minor
    title('Tx')
    xlabel('[msec]')
    subplot(2,1,2)
    t=0:Ts:Ts*(length(Signal_Rx)-1);
    plot(t/1e-3,Signal_Rx);xlabel('[msec]');grid on;grid minor
    title('Rx')
    xlabel('[msec]')
end

%% Audio

%%%%%
% [Signal_Rx,Fs]=audioread('Alon1_Rx_2.wav');
% Signal_Rx=Signal_Rx(:,2);
% load('Alon1_tx.mat')
%%%%%


%%%%%
%Signal_Tx_digital=Signal_Tx_digital*(0.95/max(abs(Signal_Tx_digital)));



% bits=16;
% recObj = audiorecorder(Fs, bits, 1); % playing object (Tx)- controls the headphones
% player = audioplayer(Signal_Tx_digital, Fs, bits); % recording object (Rx)- controls the microphone
%
% record(recObj,length(Signal_Tx_digital)*(1/Fs)*1.2); % Record your voice for T_record
% pause(0.5)
% tic
% playblocking(player); % Playing the transmitted signal
% toc
% Signal_Rx = getaudiodata(recObj);
%
%
% t_Tx=0:1/Fs:(length(Signal_Tx_digital)-1)*1/Fs;
% t_Rx=0:1/Fs:(length(Signal_Rx)-1)*1/Fs;
%
% SIGNAL_RX=fftshift(fft(Signal_Rx));
% SIGNAL_TX_DIGITAL=fftshift(fft(Signal_Tx_digital));
% figure
% subplot(3,1,1)
% plot(t_Tx,Signal_Tx_digital);grid minor
% subplot(3,1,2)
% plot(t_Rx,Signal_Rx);grid minor
% subplot(3,1,3)
% plot(db(abs(SIGNAL_RX)));grid minor
% %plot(db(abs(SIGNAL_TX_DIGITAL)));grid minor
%% Reception


t1_rx=clock;
[ Symbol_stream_Rx ] = Receiver(Signal_Rx_MIMO,Fs,OFDM_config,IF_chain_config,Simulation_config,Testing_data );
e_Rx=etime(clock,t1_rx);

%% QAM demodulation
[ data_DeCoded_Rx,data_Rx,t2_rx ] = DeModulator(Symbol_stream_Rx,Coding_config,OFDM_config,Simulation_config,Audio_config );

%% Perforamce evaluation
if strcmp(Link_Type,'HW')
    
    %Data_Rate_Audio_kbps=(Fs_Audio*q_bits)/1e3; % the audio recording produces Fs_Audio[samples/sec]*q_bits[bits/sample]  bits
    
    Data_bits=N_bits;
    e_rx=etime(t2_rx,t1_rx);
    T_communication=etime(t2_rx,t1_tx)
    Data_Rate_Comm_kbps=(Data_bits/T_communication)/1e3
    
    %     if Data_Rate_Audio_kbps>Data_Rate_Comm_kbps
    %         warning('Data rate is not sufficient')
    %     end
    
    
end

%% Analysis

hError=comm.ErrorRate;

if Coding_config.Coding_flag
    Error_strct=step(hError,data_PreCoded_Tx,data_DeCoded_Rx); % we ommit from the comparison the chopped symbols (due to group delay)
else
    Error_strct=step(hError,data_Tx,data_Rx); % we ommit from the comparison the chopped symbols: due to group delay an fgthe fact that matlab's filter function chops the end of the vector, and we connot do anything about that
end

Bits_Error=Error_strct(2);
Bits_compared_total=Error_strct(3);

if nargin<1
    
    BER=Error_strct(1)
    
    
    hEVM=comm.EVM('Normalization','Average reference signal power');
    %EVM=step(hEVM,Symbol_stream_Tx(1:length(Symbol_stream_Rx)),Symbol_stream_Rx);
    EVM=step(hEVM,Symbol_stream_Tx,Symbol_stream_Rx);
    EVM_dB=db(EVM/100)
    
    length_test=round(length(Symbol_stream_Tx)/4);
    
    EVM_start=step(hEVM,Symbol_stream_Tx(1:length_test),Symbol_stream_Rx(1:length_test));
    EVM_start_dB=db(EVM_start/100)
    
    EVM_end=step(hEVM,Symbol_stream_Tx(end-2*length_test:end-length_test),Symbol_stream_Rx(end-2*length_test:end-length_test));
    EVM_end_dB=db(EVM_end/100)
    
    %     scatterplot(Symbol_stream_Rx);
    %     grid on
    %     grid minor
    %     title(gca,['scatter plot. EVM=',num2str(EVM_dB),'[dB]'])
    %     set(gcf,'windowstyle','docked');
    
    
    Testing_data.Group_delay_Tx_total=round(Testing_data.Group_delay_total);
    
    %hCCDF = comm.CCDF('PAPROutputPort',true, 'MaximumPowerLimit', 60)  ; % adapted to the maximum range of the D/A and the A/D(-10,10)
    hCCDF = comm.CCDF('PAPROutputPort',true, 'MaximumPowerLimit', 5)  ;
   % hCCDF1=clone(hCCDF);
  %  [CCDFx,CCDFy,PAPR]=step(hCCDF,Signal_Tx_digital(Testing_data.Group_delay_Tx_total+1:end));
    [CCDFx,CCDFy,PAPR]=step(hCCDF,Signal_Tx_digital(Testing_data.Data_start_index+1:end));
    PAPR_dB=PAPR
   % PAPR1
    
    figure
    plot(hCCDF)
    title(gca,['CCDF Plot. PAPR=',num2str(PAPR_dB),'[dB]'])
    set(gcf,'windowstyle','docked')
    
    
    
    %%%H=dsp.SpectrumAnalyzer('SampleRate',Fs_air,'FrequencySpan','Span and center frequency','CenterFrequency',F_if,'Span',F_chip*4,'FrequencyResolutionMethod','RBW','RBWSource','Property','RBW',0.1,'ShowGrid',true,'YLimits',[-250,-50]);
    % H=dsp.SpectrumAnalyzer('SampleRate',Fs,'FrequencySpan','Start and stop frequencies','StartFrequency',max(0,F_if-2*F_chip),'StopFrequency',min(Fs/2,F_if+F_chip*2),'FrequencyResolutionMethod','RBW','RBWSource','Property','RBW',0.1,'ShowGrid',true,'YLimits',[-250,-50]);
    % step(H,[Signal_Tx_digital(Testing_data.Group_delay_Tx_total+1:end) Signal_Rx_digital(Testing_data.Group_delay_Tx_total+1:end)])
    % show(H)
    
    
    M=length(Signal_Rx_digital(Testing_data.Group_delay_Tx_total+1:end));
    f_adc=linspace(-Fs/2,Fs/2,M);
    f_dca=linspace(-Frec/2,Frec/2,M);
    A=fftshift(fft(Signal_Tx_digital(Testing_data.Group_delay_Tx_total+1:end),M));
    B=fftshift(fft(Signal_Rx_digital(Testing_data.Group_delay_Tx_total+1:end),M));
    
    figure
    set(gcf,'windowstyle','docked')
    plot(f_dca/1e3,db(abs(A)))
    hold on
    plot(f_adc/1e3,db(abs(B)))
    grid on
    grid minor
    xlabel('frequncy [kHz]]')
    title(['Digital IF Tx and Rx signal: at the interface of the AFE. Fchip=',num2str(F_chip/1e3),'[kHz]. Fsample=',num2str(Fs/1e3),'[kHz]. Frec=',num2str(Frec/1e3),'[kHz]'])
    legend('Tx before D/A','Rx after A/D')
    xlim([-Fs/2/1e3,Fs/2/1e3])
    
    
    Ts=1/Fs;
    Trec=1/Frec;
    t_Tx=0:Trec:Trec*(length(Signal_Tx_digital)-1);
    t_Rx=0:Ts:Ts*(length(Signal_Rx_digital)-1);
    
    P_Tx_digital=Signal_Tx_digital.^2;
    P_MA_Tx_digital=filter(ones(10000,1),1,P_Tx_digital);
    P_Rx_digital=Signal_Rx_digital.^2;
    P_MA_Rx_digital=filter(ones(10000,1),1,P_Rx_digital);
    
    figure
    set(gcf,'windowstyle','docked')
    plot(t_Tx/1e-3,db(P_MA_Tx_digital,'power'))
    hold on
    plot(t_Rx/1e-3,db(P_MA_Rx_digital,'power'))
    title('Average power at AFE s interface'); legend('Tx Signal','Rx Signal')
    grid on;grid minor
    xlabel('time [msec]');ylabel('[dB]')
    ylim([-15,-5])
    
    
    
    
    figure
    set(gcf,'windowstyle','docked')
    subplot(2,1,1)
    plot(t_Tx/1e-3,Signal_Tx_digital)
    xlabel('[msec]');grid on;grid minor
    title(['Tx Signal at AFE interface.T preamble synch=',num2str(T_preamble_synch/1e-3),'[msec]. T preamble CE=',num2str(T_preamble_CE/1e-3),'[msec].'])
    subplot(2,1,2)
    plot(t_Rx/1e-3,Signal_Rx_digital(:,1:MIMO_depth))
    xlabel('[msec]');grid on;grid minor
    title(['Rx Signal at AFE interface.T preamble synch=',num2str(T_preamble_synch/1e-3),'[msec]. T preamble CE=',num2str(T_preamble_CE/1e-3),'[msec].'])
    
    scatterplot(Symbol_stream_Rx);
    grid on
    grid minor
    title(gca,['scatter plot. EVM=',num2str(EVM_dB),'[dB]'])
    set(gcf,'windowstyle','docked');
    
    
end

end