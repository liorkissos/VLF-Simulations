function [ Signal_Out,Group_delay_total,H_channel ] = Channel(Signal_In,Delay_resp_vec,Amp_resp_vec,Fs_ref,Fs_channel,Group_delay_total)

%CHANNEL Summary of this function goes here
% Function receives a Tx real, IF carried, time domain signal, and propagates it through a
% multipath channel. AWGN is not added. based on the model
% r(t)=sum{A_n*s(t-Tau_n)} (Proakis 5th ed. equation 13.1-2)
% Channel configuration is given in OFDM modem terms, i.e; at Tchip terms.
% 
%   Detailed explanation goes here:
% Fs_ref= sampling frequency at the OFDM Modem interface, therefore,  Fs_ref= Fchip
% Fs_channel= "analog" sampling frequency. sampling frequency at the wireless channel
% Delay_resp_vec= vector of the channel's tap delays relative to the moment
% of the 1st impulse response's tap, h[0]. e.g;Delay_resp_vec=[0,4],
% Amp_resp_vec= [1, 0.5] means; h[n]=[1 0 0 0 0.5]
% Group_delay_total= sums the Tx chain's group delay with the channel
% maximal amplitude tap's delay. This item of data, the total group delay
% is necessary for a correct AWGN addition (which is a completely
% fictitious procedure, or in the case of artificial channel timing
% correction, without Minn& Zeng


%Debug_flag=1;
Debug_flag=0;


Ts_ref=1/Fs_ref;
Ts_channel=1/Fs_channel;

Delay_resp_vec=round(Delay_resp_vec*(Ts_ref/Ts_channel));

%%% Channel Creation

N=max(Delay_resp_vec+1); %Channel impulse response length- the location of the last channel's tap or the filter's length (the maximal)
n_resp_vec=Delay_resp_vec+1; % the non zero elements within the channel are in fact the delay vector elements versus the first tap plus 1. the "+1" is needed to trandform from "delay" to actual elements 
h=zeros(1,N);
h(n_resp_vec)=Amp_resp_vec;
H_channel=dfilt.df1(h,1);

%%% Signal propagation
%tic
Signal_Out=filter(H_channel,Signal_In);
%toc
if Debug_flag
    
   % fvtool(H_channel)
    
    f=linspace(-Fs_channel/2,Fs_channel/2,length(Signal_In));
    A=fftshift(fft(Signal_In));
    B=fftshift(fft(Signal_Out,length(Signal_In)));
    
    figure
    set(gcf,'windowstyle','docked')
    subplot(2,1,1)
    plot(f/1e3,db(abs(A)),f/1e3,db(abs(B)))
    grid on
    grid minor
    title('Signal spectrum before and after propgating through channel- Zoom out')
    legend(['Before'],['After'])
    subplot(2,1,2)
    plot(f/1e3,db(abs(A)),f/1e3,db(abs(B)))
    xlim([-50, 50])
    grid on
    grid minor
    title('Signal spectrum before and after propgating through channel- Zoom in')
    legend(['Before'],['After'])
   
    
    
    
    
end
%%% Group Delay

Group_delay=find(abs(h)==max(abs(h)))-1; % the channel's strongest tap delay in Ts_channel terms, relative to 1st tap (the "-1")
Group_delay_total=Group_delay_total+Group_delay;
%Group_delay_total1=Group_delay_total1+min(Delay_resp_vec);% the first nonzero element of the channel
end

