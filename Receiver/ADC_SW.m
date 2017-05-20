function [ Signal_Rx ,Fs,Testing_data] = ADC_SW( Signal_Rx,Fs,IF_chain_config,Simulation_config,Testing_data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Debug_flag=0;
%Debug_flag=1;


Group_delay_total=Testing_data.Group_delay_total;

Configuration=Simulation_config.Configuration;


%% 1) A/D & Anti-Aliasing Filter

%%% A/D: 1) Anti Aliasing Filtering 2) Sampling



%%% Anti aliasing Filter

if strcmp(Configuration,'Calibration')==0
    
   % display('the resample function limits EVM performance in high SNR conditions: ~58dB instead of ~70dB')
    
    Fs_req=IF_chain_config.Fs_req;
    
    [p,q]=rat(Fs_req/Fs);
    
    [Signal_Rx,b]=resample(Signal_Rx,p,q,40); % original value
   % [Signal_Rx,b]=resample(Signal_Rx,p,q,60);
    
    Group_delay_anti_aliasing=0; % see Help on resample: it compensates for the delay introduced
    
     Group_delay_total=Group_delay_total+Group_delay_anti_aliasing;

else
  
    N_upsample_ZOH=IF_chain_config.N_upsample_ZOH;
    
    q=N_upsample_ZOH;
    p=1;
    
    load('Filters.mat')
    
    Signal_Rx=filter(Hanti_aliasing_flt,Signal_Rx); %1st relevant sample=6645
    
    Group_delay_anti_aliasing=Hanti_aliasing_flt.order/2;
    Group_delay_total=Group_delay_total+Group_delay_anti_aliasing;
    
    d1=mod(Group_delay_total+1,q)-1;
    
    Signal_Rx=downsample(Signal_Rx,q,d1); %504
    
    Group_delay_total=Group_delay_total+mod(d1,q);

    
end

Fs=Fs*(p/q); 
Group_delay_total=Group_delay_total*(p/q); % by shifting by d1 we effectively reduce the group delay after decimation by 1, or post decimation by N_upsample_ZOH



if Debug_flag
    
    Signal_Rx_sampled=Signal_Rx;
    % we expect to see in blue N_upsample_ZOH-1 replicas of the Tx signal,
    % and in brown a sinc-shaped version of the blue curve
    
    figure(8)
    set(gcf,'windowstyle','docked');
    f=linspace(-Fs/2,Fs/2,length(Signal_Rx));
    A=fftshift(fft(Signal_Rx));
    plot(f/1e3,db(abs(A)))
    title(['Signal Tx after A/D'])
    %  legend('Signal Tx upsampled','Signal Tx upsampled & ZOH shaped','Signal Tx upsampled& ZOH shaped & reconstructed')
    xlabel('freq [kHz]')
    grid on
    grid minor
    ylim([-80 80])
    
    
%         A=Signal_Tx_upconverted(Group_delay_interp+1:length(Signal_Tx_upconverted));
%         % A=Signal_Rx_anti_aliased(Group_delay_post_Anti_Alias+1:N_upsample_ZOH:end);
%         B=Signal_Rx_sampled(Group_delay_total+1:1:end);
%     
%         EVM_Lior_Sampling=db(mean(abs(A(1:length(B))-B)/std(A)))
    
    
    
end


%% update of structures

Testing_data.Group_delay_total=Group_delay_total;

end

