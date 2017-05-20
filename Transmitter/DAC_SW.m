
function [Signal_Tx,Frec,OFDM_config,Testing_data]=DAC_SW(Signal_Tx,Frec,IF_chain_config,OFDM_config,Simulation_config,Testing_data)
%function [Signal_Tx,Frec,Group_delay_total]=DAC_SW(Signal_Tx,Frec,IF_chain_config,OFDM_config,Configuration,Group_delay_total)



Debug_flag=0;
%Debug_flag=1;

F_if=IF_chain_config.F_if;
N_upsample_ZOH=IF_chain_config.N_upsample_ZOH;

B_signal_Hz=OFDM_config.F_chip;

Group_delay_total=Testing_data.Group_delay_total;

Configuration=Simulation_config.Configuration;

%% 12) D/A & Reconstruction Filter

%%% D/A: 1) upsampling 2) filtering through ZOH filter

if (mod(N_upsample_ZOH,2)~=1 || N_upsample_ZOH<9)
    error('ZOH must have an even order and bigger than 9') % odd length to add integer number of samples to the delay. longer than 9 so that the frequency domain sinc (time domain rect) is not aliased
end

if strcmp(Configuration,'Calibration')==0 
    
    H_ZOH=ones(N_upsample_ZOH,1);
    H_ZOH_flt=dfilt.dffir(H_ZOH);
    
else
    
    load('Filters.mat')
    
end

Signal_Tx=upsample(Signal_Tx,N_upsample_ZOH); % from 491 to 5391=491*N_upsample_ZOH-(N_upsample_ZOH-1)
Signal_Tx_upsampled=Signal_Tx; % for debugging purposes

Signal_Tx=filter(H_ZOH_flt,Signal_Tx);
Signal_Tx_ZOH=Signal_Tx;

Frec=Frec*N_upsample_ZOH;
Group_delay_total=Group_delay_total*N_upsample_ZOH;



%%%% Reconstruction
if strcmp(Configuration,'Calibration')==0 % in any case other than Calibration, the interpolation factor and thus filter need to be flexible in order to adapt to the variable F_chip and the maximum Fsample cnstrained by the NI 6212
    
    Fpass=F_if+B_signal_Hz/2;
    Fstop=Frec/N_upsample_ZOH-F_if-B_signal_Hz/2; % the 1st replica's center frequency is at Fs-F_if
    Hrec_SPEC=fdesign.lowpass('N,Fp,Fst,Ap',92,Fpass*1.1,Fstop,0.001,Frec);
    Hrec_flt=design(Hrec_SPEC,'equiripple');
    
end

Signal_Tx=filter(Hrec_flt,Signal_Tx); % 1st relevant sample:6542

Group_delay_rec=Hrec_flt.order/2+((N_upsample_ZOH-1)/2); % accumulated group delay=5391+51=5442. the ZOH operation is a FILTERING (!!) by a time domain Rectangular of . Hence, we need to add to the Filter's group delay the order of the rectangular filter. emprically, need to correct the delay by -1 to (N_upsample_ZOH-1)/2-1
Group_delay_total=Group_delay_total+Group_delay_rec;
Group_delay_Tx_total=Group_delay_total; % for debugging use

%%% Power ratios
E12_1=sum(sum((abs(Signal_Tx_ZOH)).^2));
E12_2=sum(sum((abs(Signal_Tx)).^2)); % Total OFDM block (with pilots and guard bands, but w/o CP) Energy, without the preamble_CE
R12=E12_1/E12_2; %

if Debug_flag
    Signal_Tx_reconstructed=Signal_Tx;
    Group_delay_post_reconstruction=Group_delay_total;
    % we expect to see in blue N_upsample_ZOH-1 replicas of the Tx signal,
    % and in brown a sinc-shaped version of the blue curve
    figure(5)
    set(gcf,'windowstyle','docked');
    f=linspace(-Frec/2,Frec/2,length(Signal_Tx));
    A=fftshift(fft(Signal_Tx_upsampled));
    B=fftshift(fft(Signal_Tx_ZOH));
    C=fftshift(fft(Signal_Tx));
    plot(f/1e3,db(abs(A)),f/1e3,db(abs(B)),f/1e3,db(abs(C)))
    title(['Signal Tx after ZOH D/A'])
    legend('Signal Tx upsampled','Signal Tx upsampled & ZOH shaped','Signal Tx upsampled& ZOH shaped & reconstructed')
    xlabel('freq [kHz]')
    grid on
    grid minor
    ylim([-80 80])
    
    
%     A=Signal_Tx_upconverted(Group_delay_interp+1:end);
%     B=Signal_Tx_reconstructed(Group_delay_post_reconstruction+1:N_upsample_ZOH:end);
%     
%     EVM_Lior_Reconstruction=db(mean(abs(A(1:length(B))-B)/std(A)))
end


%% update of structures
OFDM_config.R12=R12;
Testing_data.Group_delay_total=Group_delay_total;

end