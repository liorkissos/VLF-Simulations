function [ Signal_Rx,SNR_other,Testing_data ] = Channel_SW( Signal_Tx,Fs,Channel_config,Simulation_config,OFDM_config,Testing_data )


%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Debug_flag=0;
%Debug_flag=1;

M=OFDM_config.M;
N_data=OFDM_config.N_data;
N_FFT=OFDM_config.N_FFT;
N_CP=OFDM_config.N_CP;

F_chip=OFDM_config.F_chip;
R0=OFDM_config.R0;
R1=OFDM_config.R1;
R2=OFDM_config.R2;
R3=OFDM_config.R3;
R4=OFDM_config.R4;
R5=OFDM_config.R5;
R6=OFDM_config.R6;
R7=OFDM_config.R7;
R10=OFDM_config.R10;
R11=OFDM_config.R11;
R12=OFDM_config.R12;
R_lengths=OFDM_config.R_lengths;


SNR_method=Channel_config.SNR_method;
SNR=Channel_config.SNR;
Delay_resp_vec=Channel_config.Delay_resp_vec; % 0=the time delay at which the first sample enters the channel
Amp_resp_vec=Channel_config.Amp_resp_vec;
Phase_resp_vec=Channel_config.Phase_resp_vec;
Fs_channel=Channel_config.Fs_channel; % the sampling frequency/interval at which the channel was given (the interval between taps)

BB_IF_flag=Simulation_config.BB_IF_flag;


Group_delay_total=Testing_data.Group_delay_total;

%% Channel %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

B_signal_Hz=F_chip; % the BB signal bandwidth. BB signal: right after CP insertion
Fs_air=Fs; % since we do not upsample yet, the air interface Fsample equals the the BB signal's bandwidth


Signal_Tx_flat=Signal_Tx;

%[ Signal_Tx,Group_delay_total,H_channel_flt ] = Channel(Signal_Tx,Delay_resp_vec,Amp_resp_vec,F_chip,Fs_air,Group_delay_total);
[ Signal_Tx,Group_delay_total,H_channel_flt ] = Channel(Signal_Tx,Delay_resp_vec,Amp_resp_vec,Fs_channel,Fs_air,Group_delay_total);


L_channel=find(H_channel_flt.numerator,1,'last' )-find(H_channel_flt.numerator, 1,'first'); % the memory (=support) of the channel is not its length but the number of taps between its first nonzero tap and its last nonzero tap, since the first zero taps will be compensated by the coarse timing synchronization

% power ratios: i case of flat channel R8=1
E8=sum((abs(Signal_Tx_flat)).^2); % Only payload subcarriers Energy. rounding is donw when group delay is not an integer
E9=sum((abs(Signal_Tx)).^2); % Total OFDM block (w/o CP) Energy
R8=E8/E9;

%%% Compensation for the group delay samples that do not contain energy:
%%% since power measurement carried out by awgn function is done over the
%%% whole signal length, including the group delay, those group delay
%%% samples may lower the average energy, thus dragging a noise addition
%%% lower than required to the data-carrying (!) samples => increase of the
%%% SNR of the data-carrying samples; EsNo_data
%%% Hence, EsNo (the actual amount sent to the awgn function)deduced from EsNo_data needs to be corrected downwards
%%% (decreased)
%%% Group delay needs to be calculated rather than evaluated in order to
%%% insert the right amount of noise


E9=sum((abs(Signal_Tx(round(Group_delay_total)+1:end))).^2); % Only payload subcarriers Energy. rounding is donw when group delay is not an integer
E10=sum((abs(Signal_Tx)).^2); % Total OFDM block (w/o CP) Energy
R9=E9/E10; % frequency domain Enery ratio (<1)

N_signal_air=length(Signal_Tx);
N_signal_air_effective=N_signal_air-Group_delay_total;


%%% SNR Addition:
%%% Important: what truly matters is whether (on calibration mode): 1)
%%% on EbNo mode: EbNo~EVM-log10(log2(M)) (e.g; M=4&EbNo=30, EVM=-33 on every signal length)
%%% 2) on EsNo_data mode: EsNo~EVM(e.g; EsNo_data=30, EVM=-30 on every constellation and every signal length)
%%% 3) BER curve is well aligned with theory

%R7=1; % fictitious; no R7 is passed from caller

switch SNR_method
    case 'EbNo' % SNR per bit on data  bits only! same as Esno_data apart from the log10(log2(M)) factor
        EbNo_data=SNR;
        EsNo_data=EbNo_data+10*log10(log2(M)); % should be equal to the resulting EVM
        EsNo=EsNo_data+10*log10(1/(R0*R1*R2*R3*R4*R5*R6*R7*R8*R10*R11*R12))+10*log10(R_lengths)+10*log10(1/R9);%10*log10((1/R8)*(N_signal_air/N_signal_air_effective));
        SNR_other=EsNo;
    case 'EsNo_data' % Calibration parameter. on Calibration mode, EVM~EsNo_data. identical to EbNo apart from being lower by log10(log2(M)) regarding EbNo for the same configuration. e.g (calibration mode); EsNo_data=30-> EVM~-30. EbNo=30-> EVM~-33
        EsNo_data=SNR;
      %  EsNo=EsNo_data+10*log10(1/(R0*R1*R2*R3*R4*R5*R6*R7*R8*R10))+10*log10(R_lengths)+10*log10(1/R9);
        EsNo=EsNo_data+10*log10(1/(R0*R1*R2*R3*R4*R5*R6*R7*R8*R10*R11*R12))+10*log10(R_lengths)+10*log10(1/R9);
        EbNo_data=EsNo_data-10*log10(log2(M)); % should be equal to the resulting EVM
        SNR_other=[]; %useful only on single link simulation. irrelevant on BER curve simuation
        
    case 'EsNo' % True Baseband signal SNR within signal's
                         %bandwidth: the SNR we would measure with a spectrum analyzer. i.e;
                        %before upsampling or IF upconversion (see Channel_SW function). EVM
                         %will eventually be different from EsNo value and
                         %will depend on length of signal.
        EsNo=SNR;
        EsNo_data=EsNo-10*log10(1/(R0*R1*R2*R3*R4*R5*R6*R7*R8*R10*R11*R12))-10*log10(R_lengths)-10*log10((1/R9));
       % EsNo_data=EsNo-10*log10(1/(R0*R1*R2*R3*R4*R5*R6*R7*R8*R10))-10*log10(R_lengths)-10*log10((1/R9));
        EbNo_data=EsNo_data-10*log10(log2(M)); % should be equal to the resulting EVM
        SNR_other=EbNo_data;
        
        
end

SNR_BB=EsNo-10*log10(Fs_air/B_signal_Hz); % takes into account the un-used samples in the time domain due to upsampling
SNR_IF=SNR_BB+3;

switch BB_IF_flag
    
    case 'IF'
        SNR_added=SNR_IF;
        
    case 'BB'
        SNR_added=SNR_BB;
end


Signal_Rx=awgn(Signal_Tx,SNR_added,'measured'); % qst relevant sample=5442


if Debug_flag
    
    fvtool(H_channel_flt)
    set(gcf,'windowstyle','docked');
    
    Signal_Rx_channel=Signal_Rx;
    % Pilots might spread along several bins becuase the bins of FFT applied over the whole
    % signal are narrower than those of the OFDM FFT, and therefor hte energy is spread along several bins
    figure(6)
    set(gcf,'windowstyle','docked');
    f=linspace(-Fs/2,Fs/2,length(Signal_Tx));
    A=fftshift(fft(Signal_Tx));
    B=fftshift(fft(Signal_Rx));
    plot(f/1e3,db(abs(A)),f/1e3,db(abs(B)))
    title(['Signal Tx after going through Channel'])
    legend('Signal Tx','Signal Rx')
    xlabel('freq [kHz]')
    grid on
    grid minor
    
    
    %     A=Signal_Tx_reconstructed(Group_delay_total+1:N_upsample_ZOH:end);
    %     B=Signal_Rx_channel(Group_delay_total+1:N_upsample_ZOH:end);
    %
    %     EVM_Lior_Channel=db(mean(abs(A(1:length(B))-B)/std(A)))
    %     display('EVM_Lior_Channel:should be tested on high SNR only')
    
end

%% update of structures

Testing_data.Group_delay_total=Group_delay_total;


end

