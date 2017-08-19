function [ Symbol_stream_Rx ] = Receiver(Signal_Rx_MIMO,Fs_ADC,OFDM_config,IF_chain_config,Simulation_config,Testing_data )
%function [ Symbol_stream_Rx,N_frames_Rx ] = Receiver(Signal_Rx,Fs,OFDM_Ref_symbol_Tx_f,Pilots_matrix_Tx,OFDM_config,IF_chain_config,Simulation_config,Testing_data )

%RECEIVER Summary of this function goes here
%   Detailed explanation goes here

%% Config
test_synch_channel_flag=0; % test only synchronization, impulse response and channel estimation
%test_synch_channel_flag=1; % test only synchronization, impulse response and channel estimation

test_signal_processing_flag=0;
%test_signal_processing_flag=1;

%% Input data

F_chip =OFDM_config.F_chip;
N_FFT = OFDM_config.N_FFT;
N_data = OFDM_config.N_data;
Npilots   =OFDM_config.Npilots;
Amp_pilots_dB  = OFDM_config.Amp_pilots_dB;
Nguard_band_left  = OFDM_config.Nguard_band_left;
Nguard_band_right  = OFDM_config.Nguard_band_right;
N_CP = OFDM_config.N_CP;
P_total=OFDM_config.P_total;

T=OFDM_config.T;
B_signal_net_Hz=OFDM_config.B_signal_net_Hz;
N_preamble_CE=OFDM_config.N_preamble_CE;
Equalizer_type=OFDM_config.Equalizer_type;
N_preamble_synch=OFDM_config.N_preamble_synch;
Enhancement_prmbl_CE=OFDM_config.Enhancement_prmbl_CE;

F_if=IF_chain_config.F_if;
Carrier_offset=IF_chain_config.Carrier_offset;
Frequency_sync_Coarse_flag=IF_chain_config.Frequency_sync_Coarse_flag;

Input_Type=Simulation_config.Input_Type;
Time_sync_flag=Simulation_config.Time_sync_flag;
BB_IF_flag=Simulation_config.BB_IF_flag;
Configuration=Simulation_config.Configuration;
Frequency_sync_flag=Simulation_config.Frequency_sync_flag;
Ref_Cfg_flag=Simulation_config.Ref_Cfg_flag;
MIMO_depth=Simulation_config.MIMO_depth;

if ~isempty(Testing_data)
    %OFDM_matrix_Tx_f_no_CP_no_Preamble=Testing_data.OFDM_matrix_Tx_f_no_CP_no_Preamble;
    Group_delay_interp=Testing_data.Group_delay_interp;
    Group_delay_total=Testing_data.Group_delay_total;
else
    Group_delay_interp=0;
    Group_delay_total=0;
end

%% Sanity checks


if(mod(N_data+1,Npilots-1)~=0)
    error('N_data_sub_carriers should be an Integer division of Npilots-1') % to enable uniform scattering of pilot subarriers
end


if(N_CP>N_FFT)
    error('N_CP should be smaller than N_FFT')
end

if mod(log2(N_FFT),1)~=0 %
    error('N_FFT should be a power of 2')
end

%% Rx Signal Demodulation %%%


for mm=1:MIMO_depth
    
    Signal_Rx=Signal_Rx_MIMO(:,mm);
    Fs=Fs_ADC;
    
    if test_signal_processing_flag
        
        Signal_Rx_sampled=Signal_Rx;
        % we expect to see in blue N_upsample_ZOH-1 replicas of the Tx signal,
        % and in brown a sinc-shaped version of the blue curve
        
        figure(8)
        set(gcf,'windowstyle','docked');
        f=linspace(-Fs/2,Fs/2,length(Signal_Rx));
        A=fftshift(fft(Signal_Rx));
        plot(f/1e3,db(abs(A)))
        title(['Signal Rx after A/D'])
        %  legend('Signal Tx upsampled','Signal Tx upsampled & ZOH shaped','Signal Tx upsampled& ZOH shaped & reconstructed')
        xlabel('freq [kHz]')
        grid on
        grid minor
        ylim([-80 80])
        
    end
    
    %% 1) Down Conversion
    
    if strcmp(BB_IF_flag,'IF')
        
        % IQ demodulation filtering will be done together with that of the
        % decimation
        
        F_if=F_if-Carrier_offset;
        
        Ts=1/Fs;
        t=0:Ts:(length(Signal_Rx)-1)*Ts;
        
        Signal_Rx=Signal_Rx.*exp(-j*2*pi*F_if*t'); %504. the value of exp(-j*2*pi*F_if*t') at n=504 should match (complex conjugate) that of Signal_Tx right after interpolation filte at n=351
        
        
        
        %%%% Phase correction: Calibration mode
        %%% Phase difference between the 2 NCO's
        
        if not(Time_sync_flag)
            Group_delay_NCO_to_NCO=Group_delay_total-Group_delay_interp;
            Phase_NCO_Tx_Rx_Acc=(2*pi*F_if)*(Group_delay_NCO_to_NCO)*Ts; % accumulated phase by the NCO Rx  until the sample arrives from the NCO Tx
            
            Signal_Rx=Signal_Rx.*exp(+j*Phase_NCO_Tx_Rx_Acc); % 504. check: (1*2*pi*F_if*t(351)-1*2*pi*F_if*t(504))+1*Phase_NCO_Tx_Rx_Acc
        end
        
        if test_signal_processing_flag
            % Pilots might spread along several bins becuase the bins of FFT applied over the whole
            % signal are narrower than those of the OFDM FFT, and therefor hte energy is spread along several bins
            figure(9)
            set(gcf,'windowstyle','docked');
            f=linspace(-Fs/2,Fs/2,length(Signal_Rx));
            A=fftshift(fft(Signal_Rx));
            plot(f/1e3,db(abs(A)))
            title(['Signal Rx after downconverting by ',num2str(F_if/1e3),' [kHz]'])
            xlabel('freq [kHz]')
            grid on
            grid minor
            
            %%% EVM not tested here versus Interpolator filter output, since an
            %%% additional filtering is needed
            
            
        end
        
    end
    
    %% 2) Decimation
    
    %%% Decimation Filter
    
    %if strcmp(%Design_flag,'Design')
    if strcmp(Configuration,'Operational')
        
        
        N_upsample=Fs/F_chip;
        
        if  Ref_Cfg_flag% operationl and reference configuration
            load('Filters_Ref_Cfg.mat')
            
        else
            
            
            if (Fs/2>2*F_if+B_signal_net_Hz)
                F_edge_replica=2*F_if-B_signal_net_Hz;
            else
                F_edge_replica=Fs-2*F_if-B_signal_net_Hz;
            end
            
            
            Fpass=B_signal_net_Hz; % actually, it is F_chip/2.  After upsampling, the right edge of the first replica decreases by N_upsample
            Fstop=min((Fs/2)/N_upsample+min(Nguard_band_left,Nguard_band_right)/T,F_edge_replica);% for no apparent reason, a higher Fstop than needed produces beter results. so, I broadened it by min(Nguard_band_left,Nguard_band_right)/T
            % Fstop=min((Fs/2)/N_upsample,F_edge_replica);%since it is a BB signal, the 2nd replica lies at an interval of Nguard_band_left/T from Fnyquist (the current Fs/2).
            
            Fpass=1.2*Fpass; % Fpass= B_signal_net_Hz is too narrow somewhow
            
            Hdec_SPEC=fdesign.lowpass('Fp,Fst,Ap,Ast',Fpass,Fstop,0.0002,70,Fs);% in the case of Nguard bad right=10
            Hdec_flt=design(Hdec_SPEC,'kaiserwin');
            
            
        end
        
        
    else % calibration
        
        load('Filters.mat')
        N_upsample=IF_chain_config.N_upsample;
        
        
    end
    
    
    if mod(Hdec_flt.order,2)~=0 && Time_sync_flag==0 % integer delay needed only when no synchronization nor equalization is made
        error('Decimating filter group delay is not an integer number')
    end
    
    Signal_Rx_downconverted=Signal_Rx;% for testing puroses
    
    if strcmp(BB_IF_flag,'IF')
        
        % Signal_Rx=2*filter(Hdec_flt,Signal_Rx); % 1st relevant sample=. the 2 compensates for the conversion loss
        Signal_Rx=2*conv(Hdec_flt.numerator,Signal_Rx,'full'); % conv instead of filter so as not to loose the suffix samples of the signal, that might cause loosing an actual payload frame
        
    else
        %Signal_Rx=1*filter(Hdec_flt,Signal_Rx); % 1st relevant sample=. the 2 compensates for the conversion loss
        Signal_Rx=1*conv(Hdec_flt.numerator,Signal_Rx,'full'); % 1st relevant sample=. the 2 compensates for the conversion loss
        
        
    end
    
    Group_delay_dec=Hdec_flt.order/2;
    Group_delay_total=Group_delay_total+Group_delay_dec;
    
    if test_signal_processing_flag
        
        Signal_Rx_downconverted_filtered=Signal_Rx;
        
    end
    
    if Time_sync_flag
        d2=0;  % no artificial offset due to prior knowledge of delay. delay will be deduced later with application of Minn& Zeng algorithm
    else
        d2=mod(mod(Group_delay_total+1,N_upsample)-1+N_upsample,N_upsample);
        
        d2=round(d2);
    end
    
    Signal_Rx=downsample(Signal_Rx,N_upsample,d2); % 1st relevant sample=
    % 21428
    
    %Group_delay_total=Group_delay_total+mod((N_upsample-d2),N_upsample);
    Group_delay_total=Group_delay_total-mod(d2,N_upsample);
    Group_delay_total=Group_delay_total/N_upsample; % by shifting by d2 we effectively reduce the group delay after decimation by 1, or post decimation by N_upsample_ZOH
    Group_delay_total_analytic=Group_delay_total;
    
    Fs=Fs/N_upsample;
    
    if test_signal_processing_flag
        
        Signal_Rx_decimated=Signal_Rx;
        
        % Pilots might spread along several bins becuase the bins of FFT applied over the whole
        % signal are narrower than those of the OFDM FFT, and therefor hte energy is spread along several bins
        
        N1=length(Signal_Rx_sampled);
        
        A=fftshift(fft(Signal_Rx_sampled,N1));
        B=fftshift(fft(Signal_Rx_downconverted,N1));
        C=fftshift(fft(Signal_Rx_downconverted_filtered,N1));
        D=fftshift(fft(Signal_Rx_decimated,N1));
        
        
        figure(10)
        set(gcf,'windowstyle','docked');
        subplot(3,1,1)
        f=linspace(-Fs*N_upsample/2,Fs*N_upsample/2,N1);
        plot(f/1e3,db(abs(A)),f/1e3,db(abs(B)))
        title(['Signal Rx before and after Downconvertion by ',num2str(F_if/1e3),' [kHz]'])
        xlabel('freq [kHz]')
        grid on
        grid minor
        legend(['after'],['before'])
        ylim([-60 40])
        
        subplot(3,1,2)
        plot(f/1e3,db(abs(B)))
        hold on
        plot(f/1e3,db(abs(C)))
        title(['Signal Rx after Downconversion by ',num2str(F_if/1e3),'[kHz], and filtering'])
        xlabel('freq [kHz]')
        grid on
        grid minor
        legend(['after'],['before'])
        ylim([-100 50])
        
        
        subplot(3,1,3)
        f=linspace(-Fs/2,Fs/2,N1);
        plot(f/1e3,db(abs(D)))
        title(['Signal Rx after Downconversion by ',num2str(F_if/1e3),'[kHz],filtering  and Decimation by ',num2str(N_upsample),' '])
        xlabel('freq [kHz]')
        grid on
        grid minor
        ylim([-100 50])
        
        %  A=Signal_Tx_OFDM;
        %  B=Signal_Rx_decimated(Group_delay_total+1:end);
        
        %     EVM_Lior_Decimation=db(mean(abs(A(1:length(B))-B)/std(A)))
        %     EVM_Lior_Decimation1=db(mean(abs(A(33:156)-B(33:156))/std(A(33:156))))
        
    end
    
    %% OFDM: Time domain %%%
    
    %% 3) Timing & Carrier Offsets: Metrics calculation
    
    % Timing synchronization follows Minn& Zeng Algorithm :"On Timing Offset
    % Estimation for OFDM Systems". 4 estimators based on the maximization of 4 likelihood functions, are suggested there, and we chose the method B with modification #1,
    % as recommended on the last line of the article.
    % That algorithm synchronizes only at Tchip resolution, which is quite coarse.
    % The delay is likely to contain a fractional part of a sample (sampled at
    % Tchip), which effects the phase of the sample. that part will be
    % compensated by the equalizer
    
    
    m=log2(N_preamble_synch*N_FFT/4); % Minn & Zeng. L=N_FFT/4
    N_preamble_symbols=(2^m-1)/log2(2); % preamble is always BPSK, thus log2(2)
    
    L=N_preamble_symbols; % L should verify: N_FFT=4*L+4
    
    
    % Loop Initializations
    Metric_sync=zeros(1,length(Signal_Rx));
    P_time=zeros(1,length(Signal_Rx));
    R=zeros(1,length(Signal_Rx));
    P_frequency=zeros(1,length(Signal_Rx));
    
    nn=1;
    Metric_sync_max=0.85; % high enough to avoid detecting false peaks
    nn_max=length(Signal_Rx);
    d_mp=min(25,round(N_CP/2)); % the parameter that defines the clean part of the subsequnce. can be as high as N_CP, the maximum length channel presumably correctable by OFDM
    while nn+4*L-1<length(Signal_Rx) % nn+4*L-1 is the farmost index of the 4th quarter of examined portion of the signal
        
        
        % Portions of received signal;associated with the subsequences of the preamble_synch
        % generated in the transmitter. will be used for timing synchronization
        % as well as coarse frequency alignment
        P1=Signal_Rx(nn:L+nn-1);
        P2=Signal_Rx(L+nn:2*L+nn-1);
        P3=Signal_Rx(2*L+nn:3*L+nn-1);
        P4=Signal_Rx(3*L+nn:4*L+nn-1);
        
        %%%% Timing synchronization metric calculation
        
        % 1)the Nominator of the Likelihood function: equation 11 in the article
        P_time(nn)=P1'*P2+P3'*P4; % identical to P(nn)=Signal_Rx(nn:L+nn-1)'*Signal_Rx(L+nn:2*L+nn-1)+Signal_Rx(2*L+nn:3*L+nn-1)'*Signal_Rx(3*L+nn:4*L+nn-1);
        
        % 2) the Denominator of the likelihood function: equation 8 in the article
        % (that is the modification. original method suggested equation 12 as
        % the denominator
        R(nn)=(1/2)*(Signal_Rx(nn:4*L+nn-1)'*Signal_Rx(nn:4*L+nn-1));
        % 3) The Metric that has better performance among the 2 proposed (see
        % summary of article)
        Metric_sync(nn)=(P_time(nn)'*P_time(nn))/(R(nn).^2);
        
        %%% Premature peak detection: detect the required peak without going
        %%% all over the signal. since it is normalized, it should be little
        %%% less than 1. the exact value depends on the multipath pattern (the more flat the closer to 1)
        if nn>2
            
            Metric_sync_max=max(Metric_sync_max,Metric_sync(nn));
            Metric_sync_buffer=Metric_sync(nn-2:nn); % 3 terms length sub vector
            
            if sign(Metric_sync_buffer-Metric_sync_max)==[-1 0 -1] % means that peak is in the middle point
                index_Metric_max=nn-1; % the middle point
                nn_max=index_Metric_max+80; % continue simulation beyond maximum to enable detecting a possible higher peak and because P_frequency's peak usually comes later
                
            end
            
        end
        
        %%%% Frequency syncronization metric calculation; based on my own idea:
        %%%% wo multipath and noise, all 4 subsequences should be similar. in a
        %%%% multipath channel the beginning and the end of each sequence are
        %%%% contaminated and hence different from one subsequence to the
        %%%% other. the idea is to correlate only the "clean" part of the
        %%%% subsequences; its length is shorter the longer the channel's
        %%%% memory is.
        %%%% The shorter the clean parts are, the smaller the correlation is,
        %%%% thus more susceptible to be affected by noise. Therefore, we chose
        %%%% the preamble_synch to be a sequence of (N_preamble_synch*N_FFT+N_CP)*Tchip
        %%%% duration (while a standard symbol's duration is (1*N_FFT+N_CP)*Tchip
        P_frequency(nn)=P1(d_mp:end-d_mp)'*P2(d_mp:end-d_mp)+P3(d_mp:end-d_mp)'*P4(d_mp:end-d_mp);
        
        
        if nn==nn_max
            break
        else
            nn=nn+1;
        end
        
        
        %  nn=nn+1;
    end
    
    if nn==length(Signal_Rx)
        error('premature peak detection failed') % detection of the peak without running all along the signal
    end
    
    %% 4) Timing Offset correction
    
    if Time_sync_flag % timing synchronization based on the self calculation with the Minn&Zeng algorithm. else: group delay is the one calculated based on
        Group_delay_total=find(Metric_sync==max(Metric_sync))-N_CP-1;
        
        if strcmp(Configuration,'Impulse Response')
            % when running an impulse OFDM block (see transmitter part), there happens to
            %be 2 (instead of 1)occasions where we get a maxima. the 1st is the
            %true one, and the 2nd happens right at the middle of the frame,
            %because the denominator of the metric, R, has a minima. Therefore,
            %we choose the 1st, the minimum
            
            Group_delay_total=min(find(Metric_sync>0.99*max(Metric_sync)))-N_CP-1;
            
            
        end
        
        % Deliberate delay: read summary on the issue. needed in order to create a deliberate timing error in direction of the CP (backwards).
        % this will assure that 1st sample of CP is included on one hand. on
        % the other hand, that error is fully compensable as long as it is
        % lower than N_CP
        
        if N_FFT==512
            offset_deliberate=round(1*N_CP/5);
        else
            offset_deliberate=round(2*N_CP/5);
        end
        
        if offset_deliberate>N_CP-1;
            error('deliberate synchronization offset must not be greater than N_CP')
        end
        
        Group_delay_total=Group_delay_total-offset_deliberate; %% greatly improves results. why???? because of causality of digital processing?
        
    end
    
    
    
    if  test_synch_channel_flag
        %    P_abs=abs(P);
        figure(11)
        set(gcf,'windowstyle','docked');
        plot(Metric_sync)
        title('Minn& Zeng Likelihood function: method B, modification # 1 ')
        hold off
        
        figure(12)
        set(gcf,'windowstyle','docked');
        plot(abs(P_time))
        hold on
        plot(R)
        
        Group_delay_total_estimated=Group_delay_total
        Group_delay_total_analytic
        
        %%% calculation of the memory (=support) of the channel seen by the
        %%% pure OFDM Modem. i.e; everything between the S/P and the P/S
        
        %%% Individual filters support. downsampling operations are transfered
        %%% backwards, towards Tx chain so that it will be cancelled by the
        %%% upsampling operations.
        %%% length of filter (H(z))= order of filter (H(z))+1, say N
        %%% length of filter (H(z^1/M))= ceil(N/M)
        %     Linterp=ceil((Hinterp_flt.order+1)/N_upsample);
        %     L_inv_sinc=ceil((H_inv_sinc_flt.order+1)/N_upsample);
        %     L_ZOH=ceil((H_ZOH_flt.order+1)/(N_upsample*N_upsample_ZOH));
        %     Lrec=ceil((Hrec_flt.order+1)/(N_upsample*N_upsample_ZOH));
        %     L_channel=ceil(L_channel/(N_upsample*N_upsample_ZOH));
        %     Lanti_aliasing=ceil((Hanti_aliasing_flt.order+1)/(N_upsample*N_upsample_ZOH));
        %     Ldec=ceil((Hdec_flt.order+1)/N_upsample);
        %
        %     %%% Total support: the support of the concatenation of 2 filters, with
        %     %%% lengths of M and L accordingly is M+L-1
        %     L_total=Linterp+L_inv_sinc-1+L_ZOH-1+Lrec-1+L_channel-1+Lanti_aliasing-1+Ldec-1
        
        %     if L_total>N_CP
        %         warning('Channel memory is greater than N_CP')
        %     end
        
        %     A=Signal_Rx(Group_delay_total+1:end);
        %     B=A;
    end
    
    %% 4) Carrier Offset correction: Coarse # 1
    
    T_chip=1/F_chip;% Sampling time of the OFDM transmitter/Reciver, before/after the Front-End
    
    if Frequency_sync_flag && Frequency_sync_Coarse_flag % coarse carrier synch harms performance if no Carrier offset is present. especially on low SNR
        
        %%% finding the frequency offset
        
        %     %%%% TEMP
        %     find(abs(P_frequency)==max(abs(P_frequency)))
        %    find(abs(Metric_sync)==max(abs(Metric_sync)))
        %     %%%% TEMP
        
        P_frequency_max=P_frequency(abs(P_frequency)==max(abs(P_frequency)));
        Phase_offset=angle(P_frequency_max);
        Carrier_offset_est_coarse=Phase_offset/(2*pi*(L)*T_chip);
        
        %%% testing the correctness
        % Carrier_offset
        
        Carrier_offset_est_error_dB=db(abs(Carrier_offset-Carrier_offset_est_coarse)/Carrier_offset);
        
        %%% the actual frequency correction
        t=0:T_chip:(length(Signal_Rx)-1)*T_chip;
        Signal_Rx=Signal_Rx.*exp(-j*2*pi*Carrier_offset_est_coarse*t');
        
    else
        Carrier_offset_est_coarse=0;
        
        
    end
    
    %% 5) Chopping group delay
    
    Signal_Rx=Signal_Rx(Group_delay_total+1:end); %15300: cal, 14049: op
    
    %% 6) Preamble_Synch removal
    Signal_Rx=Signal_Rx(N_CP+(N_preamble_synch)*N_FFT+1:end);
    
    %% 7) Serial to Parralel
    
    %%% Re-arrangement: removal of the group delay samples and chopping the end samples to fit
    %%% to an integer number of OFDM frames
    
    if mod(Group_delay_total,1)~=0
        error('non integer group delay. need to modify filters lengths in order to create integer divisions with upsampling factors')
    end
    
    
    % Signal_Rx1=Signal_Rx;
    % % %%% OLD -start
    % N_frames_Rx=N_frames_Rx-N_preamble_synch;
    % Signal_Rx=Signal_Rx(1:(N_FFT+N_CP)*N_frames_Rx);
    % %%% OLD- end
    
    % %%% NEW -start
    tail=mod(length(Signal_Rx),N_FFT+N_CP);%
    N_frames_Rx=(length(Signal_Rx)-tail)/(N_FFT+N_CP);
    Signal_Rx=Signal_Rx(1:(N_FFT+N_CP)*N_frames_Rx);
    % %%% NEW -end
    
    
    %N_frames_Rx=floor(length(Signal_Rx)/(N_FFT+N_CP)); %
    Signal_matrix_Rx=reshape(Signal_Rx,[N_FFT+N_CP,N_frames_Rx]);
    
    %% 8) CP removal
    
    OFDM_matrix_Rx_t=Signal_matrix_Rx(N_CP+1:end,:);
    
    %% 8*) Impulse response (test mode only)
    
    if test_synch_channel_flag
        
        if strcmp(Input_Type,'Impulse') % && Time_sync_flag==0 % no use in testing the channel before
            
            OFDM_matrix_Rx_t=OFDM_matrix_Rx_t(:,2:end); %there are 2 preamble symbols in the case input=Impulse. thus another preamble has to be removed
            
            
            Channel_resp_t_1=OFDM_matrix_Rx_t(:,1);
            Channel_resp_t_2=OFDM_matrix_Rx_t(:,2);
            
            
            Channel_resp1=fftshift(fft(Channel_resp_t_1,N_data));
            Channel_resp2=fftshift(fft(Channel_resp_t_2,N_data));
            
            
            figure(13)
            set(gcf,'windowstyle','docked');
            subplot(3,1,1)
            stem(2:length(Channel_resp_t_1),real(Channel_resp_t_1(2:end)))
            hold on
            stem(2:length(Channel_resp_t_2),real(Channel_resp_t_2(2:end)))
            hold off
            title('Time Domain:Impulse response- chip number 2 (!) and on')
            xlabel('chip')
            legend('Symbol #1','Symbol #2')
            grid on
            grid minor
            xlim([1 L_total]) % L_total= memory of the effective channel
            
            subplot(3,1,2)
            plot(angle(Channel_resp1)/pi)
            hold on
            plot(angle(Channel_resp2)/pi)
            hold off
            title(['Frequency Domain: Phase. N data=',num2str(N_data),' subcarriers'])
            xlabel('subcarrier')
            ylabel('angle/pi')
            legend('Symbol #1','Symbol #2')
            grid on
            grid minor
            
            subplot(3,1,3)
            plot(db(abs(Channel_resp1)))
            hold on
            plot(db(abs(Channel_resp2)))
            hold off
            title(['Frequency Domain: Amplitude. N data=',num2str(N_data),' subcarriers'])
            xlabel('subcarrier')
            ylabel('dB')
            legend('Symbol #1','Symbol #2')
            grid on
            grid minor
            
            return
            
        end
        
    end
    
    %% 9) Preamble_CE power attenuation:
    %%% no need to attenuate synchronization preamble,only channel estimation
    %%% preamble, because not compared to anything in the receiver (unlike the
    %%% channel estimation preamble)
    
    
    
    %OFDM_matrix_Rx_t(:,1:N_preamble_CE)=OFDM_matrix_Rx_t(:,1:N_preamble_CE)*db2mag(-3.5); % attenuation of the preamble_CE power enhancement
    OFDM_matrix_Rx_t(:,1:N_preamble_CE)=OFDM_matrix_Rx_t(:,1:N_preamble_CE)*db2mag(-Enhancement_prmbl_CE); % attenuation of the preamble_CE power enhancement
    
    %% 9) DFT
    
    OFDM_matrix_Rx_f=fft(OFDM_matrix_Rx_t);
    OFDM_matrix_Rx_f=fftshift(OFDM_matrix_Rx_f,1); % fftshift over the columns. brings them to the original arrangement of [-N_FFT/2,N_FFT/2] instead of [0,N_FFT]
    
    %% OFDM: Frequency domain %%%
    
    %% 10) Guard Bands removal
    
    OFDM_matrix_Rx_f=OFDM_matrix_Rx_f(Nguard_band_left+1:end-Nguard_band_right,:); % removing only guard bands. not pilots
    
    %% 11) Equalizer
    
    switch N_FFT % sequences include the DC subcarrier
        case 512
            load('PN_seq_512.mat')
        case 256
            load('PN_seq_256.mat') %%% Sequence taken from: IEEE C802.16a-02/17 document, which contains particularly low PAPR sequences
        case 64
            %   load('PN_seq_64.mat')
            load('PN_seq_64_1.mat')
    end
    
    Preamble_CE_stream_Tx=AAA.';
    
    P_preamble_CE=P_total/((N_data)+Npilots); % the subcarriers in this symbol all hav ethe same power
    
    OFDM_Ref_symbol_Tx_f=Preamble_CE_stream_Tx*sqrt(P_preamble_CE/1); % scaling the preamble_CE subcarriers to P_preamble_CE power; the average power per premable_CE subcarrier. before scaling, since these are BPSK symbols, their power is (+/-1)^2=1
    
    if Time_sync_flag
        
        OFDM_Ref_symbol_Rx_f=OFDM_matrix_Rx_f(:,1); % preamble has been removed. 1st symbol (1st column) is no longer preable
        
        switch Equalizer_type
            
            case 'LS'
                %OFDM_Ref_symbol_Rx_f=OFDM_matrix_Rx_f(:,1); % preamble has been removed. 1st symbol (1st column) is no longer preable
                
                H_est=OFDM_Ref_symbol_Rx_f./OFDM_Ref_symbol_Tx_f;
                
                if test_synch_channel_flag
                    
                    test_symbol_1=3;
                    test_symbol_2=4;
                    
                    
                    OFDM_Ref_symbol_Tx_f1=OFDM_matrix_Tx_f_no_CP_no_Preamble(:,test_symbol_1); % Refernce symbol for the purpose of equalization is the 1st
                    OFDM_Ref_symbol_Tx_f2=OFDM_matrix_Tx_f_no_CP_no_Preamble(:,test_symbol_2); % Refernce symbol for the purpose of equalization is the 1st
                    
                    
                    OFDM_Ref_symbol_Rx_f1=OFDM_matrix_Rx_f(:,test_symbol_1); % preamble has been removed. 1st symbol (1st column) is no longer preable
                    OFDM_Ref_symbol_Rx_f2=OFDM_matrix_Rx_f(:,test_symbol_2); % preamble has been removed. 1st symbol (1st column) is no longer preable
                    
                    
                    H_est1=OFDM_Ref_symbol_Rx_f1./OFDM_Ref_symbol_Tx_f1;
                    H_est2=OFDM_Ref_symbol_Rx_f2./OFDM_Ref_symbol_Tx_f2;
                    
                    figure(13)
                    set(gcf,'windowstyle','docked')
                    subplot(2,1,1)
                    plot(angle(H_est)/pi)
                    hold on; plot(angle(H_est1)/pi)
                    hold on; plot(angle(H_est2)/pi)
                    title(['equalizer angle. N data=',num2str(N_data),' subcarriers'])
                    xlabel('subcarrier')
                    ylabel('angle/pi')
                    grid on
                    grid minor
                    %legend('upon 1st symbol','upon 2nd symbol','upon 3rd symbol')
                    legend(['upon 1st symbol'],['upon ',num2str(test_symbol_1),'th symbol'],['upon ',num2str(test_symbol_2),'th symbol'])
                    
                    
                    subplot(2,1,2)
                    plot(abs(H_est))
                    hold on; plot(abs(H_est1))
                    hold on; plot(abs(H_est2))
                    title(['equalizer norm. N data=',num2str(N_data),' subcarriers'])
                    xlabel('subcarrier')
                    %ylabel('')
                    grid on
                    grid minor
                    legend(['upon 1st symbol'],['upon ',num2str(test_symbol_1),'th symbol'],['upon ',num2str(test_symbol_2),'th symbol'])
                end
                
                % H_est_mat=repmat(H_est,[1,N_frames_Rx]);
                
                % h_est=inv(B)*H_est; % do not erase! needs to be debugged.
                % should realize an IDFT
                
            case 'ML'
                
                L=N_CP;
                D=N_FFT-(Nguard_band_left+Nguard_band_right+1);
                
                W_twiddle=@(km,l,L) (exp(-j*(2*pi/N_FFT)*km'*l)); % BFB, p.718, equation (20.13)
                km_vec=[-N_FFT/2+Nguard_band_left:-1,1:(N_FFT/2-1)-Nguard_band_right]; % only data carrying subcarriers. we skip the DC subcarrier, index=0
                l_vec=0:1:L-1;
                
                B=W_twiddle(km_vec,l_vec,L);% performs the operation: T1=exp(-j*(2*pi/L)*(0:L-1)'*(0:L-1));
                
                S=diag([OFDM_Ref_symbol_Tx_f(1:D/2);OFDM_Ref_symbol_Tx_f(D/2+2:end)]); % BFB, p 718 between 20.13 and 20.14
                z=[OFDM_Ref_symbol_Rx_f(1:D/2);OFDM_Ref_symbol_Rx_f(D/2+2:end)]; % Symbol without Guard bands. (1:26) and (28:53) in case of N_FFT=64
                
                h_est=(inv(B'*(S')*S*B))*(B'*S'*z); % BFB, (20.16)
                H_est=B*h_est; % result is a D-1 long vector. B applies the DFT transform onto h. BFB (20.13) page 718.
                
                
                
                if test_signal_processing_flag
                    figure(1)
                    set(gcf,'windowstyle','docked')
                    stem(real(h_est))
                    hold on
                    stem(imag(h_est))
                    title('Channe estimation- Time domain')
                    legend(['real'],['imaginary'])
                    
                    figure(2)
                    set(gcf,'windowstyle','docked')
                    
                    %     stem(db(abs(H_est1)))
                    %     hold on
                    stem(db(abs(H_est)))
                    title('Channe estimation- Frequency domain')
                    
                end
                
                dummy=0.1; % a fictitious DC subcarrier. only for the sake of being compatible with LS equalizer which was designed previously
                H_est=[H_est(1:D/2);dummy;H_est(D/2+1:D)];
              %  H_est_mat=repmat(H_est,[1,N_frames_Rx]); %%% TEMP 24/2
                
                
        end
        
       %  OFDM_matrix_Rx_f=OFDM_matrix_Rx_f./H_est_mat; % the equalization operation  %%% TEMP 24/2
       OFDM_matrix_Rx_f_MIMO(:,:,mm)=OFDM_matrix_Rx_f;
       H_est_MIMO(:,mm)=H_est;
    end
    
end % MIMO loop

%% MIMO: Maximum Ratio Combination


if Time_sync_flag
    
    H_est_MIMO_norm=diag(H_est_MIMO*H_est_MIMO');
    W_estimator_MIMO=H_est_MIMO./repmat(H_est_MIMO_norm,[1,MIMO_depth]); % BFB, p.60, 735 (eq. 735)

    for kk=1:N_data+1+Npilots
        OFDM_subcarrier_kk_Rx_MIMO=squeeze(OFDM_matrix_Rx_f_MIMO(kk,:,:));
        W_estimator_subcarrier_kk_MIMO=W_estimator_MIMO(kk,:);
        
        OFDM_matrix_Rx_f(kk,:)=conj(W_estimator_subcarrier_kk_MIMO)*OFDM_subcarrier_kk_Rx_MIMO.'; % BFB p.734, eq.20.75

    end

end

%% 12) Noise variance estimation

u_vec=OFDM_matrix_Rx_f(:,1)-OFDM_matrix_Rx_f(:,2); % BFB, p.721, equation (20.28). supposing they are both Preamble_CE

Ruu=u_vec*u_vec';

Rnoise=Ruu/2; % BFB, p.721, (20.30)

%% 13) Carrier offset synchronization- Coarse # 2 (not realized yet)

if 0
    if Frequency_sync_flag && Frequency_sync_Coarse_flag % coarse carrier synch harms performance if no Carrier offset is present. especially on low SNR
        
        Phase_offset_est_coarse_vec=angle(OFDM_matrix_Rx_f(:,1))-angle(OFDM_Ref_symbol_Tx_f);
        
        Phase_offset_est_coarse_vec=mod(Phase_offset_est_coarse_vec+2*pi,2*pi);
        
        Phase_offset_est_coarse=mean(Phase_offset_est_coarse_vec,1);
        
        % extention: not really needed- equivalent frequency error causing the
        % phase shifts above. no actual use
        Carrier_offset_est_coarse=Phase_offset_est_coarse./(2*pi*T_chip*((0)*(N_FFT+N_CP)+0.5*N_FFT)); % Prasad, p.131, 5.23: the phase part
        Carrier_offset_est=Carrier_offset_est_coarse-mean(Carrier_offset_est_coarse(2:end));
        Carrier_offset_est_error_dB=db(abs(Carrier_offset_est-Carrier_offset)/Carrier_offset);
        
        %%% the calculation's result is not being used to correct anything
        
    end
end

%% 14) Preamble CE removal

%N_preamble_CE=10;
% removal of the channel estimtion preamble
OFDM_matrix_Rx_f=OFDM_matrix_Rx_f(:,N_preamble_CE+1:end); % removal of the channel estimtion preamble
N_frames_Rx=N_frames_Rx-N_preamble_CE; % reflects the removal of the 2 preamble_CE's (=channel estimating)


%% IPTS

if OFDM_config.PTS.PTS
    
    [ OFDM_matrix_Rx_f ] = IPTS_Rx( OFDM_matrix_Rx_f,OFDM_config);
    
end

%% 15) Pilots removal

if Npilots % if pilots are present
    
    chunk_size=(N_data+1)/(Npilots-1);
    Pilots_matrix_Rx(1,:)=OFDM_matrix_Rx_f(1,:);
    
    OFDM_matrix_Rx_f=OFDM_matrix_Rx_f(2:end,:); % removal of the leftmost pilot (currently- 1st row of matrix)
    
    nn=1;
    data_subcarrier_matrix_Rx=[];
    for kk=1:Npilots-1
        chunk=OFDM_matrix_Rx_f(nn:nn+chunk_size-1,:);
        data_subcarrier_matrix_Rx=cat(1,data_subcarrier_matrix_Rx,chunk);
        
        Pilots_matrix_Rx(kk+1,:)=OFDM_matrix_Rx_f(nn+chunk_size,:);
        
        nn=nn+chunk_size+1;
    end
    
    place=(N_FFT/2+1)-(Nguard_band_left+Npilots/2)-1; % the place of the DC carrier
    data_subcarrier_matrix_Rx(place+1,:)=[];
    
else
    
    data_subcarrier_matrix_Rx=OFDM_matrix_Rx_f;
    
end

% if save_SC_matrix
%     data_subcarrier_matrix_Tx_Dan=data_subcarrier_matrix_Tx(:,1:N_frames_Rx);
%     save(['data subcarrier matrix Tx. SNR=',num2str(SNR),'.mat'], 'data_subcarrier_matrix_Tx_Dan' )
%
%     data_subcarrier_matrix_Rx_Dan=data_subcarrier_matrix_Rx;
%     save(['data subcarrier matrix Rx. SNR=',num2str(SNR),'.mat'], 'data_subcarrier_matrix_Rx_Dan' )
%
% end

%% 16) Carrier offset synchronization- Fine

%%% In fact, this block is realized more as a Phase synchronizer than as
%%% fine frequency synchronizer


Amp_pilots=db2pow(Amp_pilots_dB);
P_data=P_total/(Amp_pilots*Npilots+N_data); % average power of data subcarrier


Pilots_matrix_Tx=sqrt(Amp_pilots*P_data)*ones(Npilots,N_frames_Rx);



if Frequency_sync_flag % always operational
    
    %Phase_offset_est_fine=angle(Pilots_matrix_Tx(:,1:N_frames_Rx))-angle(Pilots_matrix_Rx);
    Phase_offset_est_fine=angle(Pilots_matrix_Tx)-angle(Pilots_matrix_Rx);
    
    Phase_offset_est_fine=mean(Phase_offset_est_fine,1);
    
    Phase_offset_est_fine_mat=repmat(Phase_offset_est_fine,[N_data,1]);
    
    
    data_subcarrier_matrix_Rx=data_subcarrier_matrix_Rx.*exp(j*Phase_offset_est_fine_mat);% the actual correction of the phase, symbol (column) by symbol
    
    % extention: not really needed- equivalent frequency error causing the
    % phase shifts above. no actual use
    Carrier_offset_est_fine=Phase_offset_est_fine./(2*pi*T_chip*((0:N_frames_Rx-1)*(N_FFT+N_CP)+0.5*N_FFT)); % Prasad, p.131, 5.23: the phase part
    Carrier_offset_est=Carrier_offset_est_coarse-mean(Carrier_offset_est_fine(2:end));
    Carrier_offset_est_error_dB=db(abs(Carrier_offset_est-Carrier_offset)/Carrier_offset);
    
    
    
end

%%% Estimation of remainder of frequency error: Prasad p.131. equation 5.23:
%%% through solving a LS problem: phi_vec=A_mat*x_vec. remaining frequency
%%% offset, additionally to harming phase, harms amplitude too by shaping
%%% it in form of sinc. the problem is more severe if frequency offset
%%% varies with time, making the initial estimation invalid.
%%% where: phi(k)=teta+delta_f*(2*pi*k*T+pi*T_FFT). x_vec=[teta, delta_f]
if Frequency_sync_flag && test_synch_channel_flag % Prasad p.131. equation 5.23: estimation of remainder of frequency error
    
    MMM=16; % number of splits of the frame; controls the over-determination of the LS problem. the longer the splits (the smaller the MMM), the more overdetermoned the problem is and the better the estmation is
    K=floor(N_frames_Rx/MMM); % the sub-frame's length on super symbol terms
    
    if K<45
        warning('estimation is not sufficiently over determined. MMM is too large')
    end
    
    delta_f_old=0;
    format long
    for kk=1:MMM
        phi_vec=Phase_offset_est_fine(1+(kk-1)*K:min(N_frames_Rx,kk*K))'; % empiric phase offsets:  measured previously uponpilot subcarriers
        A_mat=[ones(K,1),2*pi*(N_FFT+N_CP)*T_chip*(1:K)'+pi*N_FFT*T_chip]; % the LS problem matrix
        x_vec=inv(A_mat'*A_mat)*A_mat'*phi_vec; % x_vec(1)=teta, constant phase. x_vec(2)=remaining frequency offset (after time domain correction)
        x_vec_disp(:,kk)=x_vec;
        
        
        %%%
        
        delta_f=x_vec(2);
        delta_f_avg=(delta_f_old+delta_f)/2;
        correction1=1/sinc(delta_f_avg*N_FFT*T_chip)
        correction=1/diric(delta_f_avg*N_FFT*T_chip,8)
        data_subcarrier_matrix_Rx(1+(kk-1)*K:min(N_frames_Rx,kk*K))=correction*data_subcarrier_matrix_Rx(1+(kk-1)*K:min(N_frames_Rx,kk*K));
        
        delta_f_old=delta_f;
        
        %%%
        
        %%%% if we ever want to correc the phase this way, uncomment next 3
        %%%% lines (same results as simple phase correction already
        %%%% implemeted)
        %         Phase_offset_est_fine1=[x_vec(1)*ones(K,1)+x_vec(2)*(2*pi*(1:K)'*(N_CP+N_FFT)*T_chip+pi*N_FFT*T_chip)];
        %         Phase_offset_est_fine_mat=repmat(Phase_offset_est_fine1',[N_data,1]);
        %         data_subcarrier_matrix_Rx(:,1+(kk-1)*K:min(end,kk*K))=data_subcarrier_matrix_Rx(:,1+(kk-1)*K:min(end,kk*K)).*exp(j*Phase_offset_est_fine_mat);
    end
    format short
    
    figure
    plot(x_vec_disp(2,:))
    set(gcf,'windowstyle','docked')
    grid minor
    title('Carrier offset Vs. the super symbols of the frame')
    ylabel('Carrier offset')
    xlabel('Super Symbol index within the frame')
    
    
    
    
end

%% 17) QAM symbols Matrix to QAM symbols Vector

Symbol_stream_Rx=reshape(data_subcarrier_matrix_Rx,[(N_data)*N_frames_Rx,1]);

end

