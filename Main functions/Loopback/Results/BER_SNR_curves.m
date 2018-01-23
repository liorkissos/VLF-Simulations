
clc
clear
close all

% %%% Waveform choosing
% 
% M=16;
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.95,2.95. Coding=1. Code rate=0.73333.mat')
% 
% addition=2;
% EbNo_data_vec=[EbNo_data_vec,max(EbNo_data_vec)+addition];
% BER_Ref = berawgn(EbNo_data_vec,'qam',max(4,M));
% EsNo_vec=[EsNo_vec,EsNo_vec(end)+addition];
% 
% 
% 
% h_figure=figure(1);
% set(h_figure,'WindowStyle','Docked')
% 
% semilogy(EsNo_vec,BER_Ref,'-o')
% hold on
% grid on
% grid minor
% xlabel('EsNo [dB]')
% ylabel('BER')
% xlim([0,max(EsNo_vec)])
% ylim(gca,[0.99e-5 1])
% %legend('')
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.95,2.95. Coding=1. Code rate=0.73333.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.95,2.95. Coding=1. Code rate=0.6.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.95,2.95. Coding=1. Code rate=0.46667.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% 
% legend('no coding','coding rate=0.73','coding rate=0.6','coding rate=0.48')
% 

%%%  Metric Synch

% Fchip=10e3;
% Tchip=1/Fchip;
% 
% 
% 
% load('High SNR, No  Multipath.mat')
% 
% t=0:Tchip:(length(Metric_sync)-1)*Tchip;


% figure(11)
% %set(gcf,'windowstyle','docked');
% plot(t/Tchip,Metric_sync)
% grid minor
% title(['Minn& Zeng Likelihood function: method B, modification # 1. Fchip=10 kHz '])
% xlabel('Time [Tchip]')
% ylim ([0 1.1])
% 
% hold on
% load('Low SNR, No  Multipath.mat')
% plot(t/Tchip,Metric_sync)
% 
% hold on
% load('Low SNR, With  Multipath.mat')
% plot(t/Tchip,Metric_sync)
% 
% hold on
% load('Low SNR, with  Multipath, short preamble.mat')
% plot(Metric_sync)
% 
% legend('High SNR, No  Multipath','Low SNR, No  Multipath','Low SNR, With  Multipath','Low SNR, With  Multipath, Short Preamble')
% 

% %%% Equalizer type
% 
% M=16;
% %load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=LSChannel=MP 2 taps 18 chip 0.5.mat')
%  load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.95,2.95. Coding=1. Code rate=0.73333.mat')
% 
% 
% 
%  addition=2;
%  EbNo_data_vec=[EbNo_data_vec,EbNo_data_vec(end)+2,EbNo_data_vec(end)+3];
%  BER_Ref = berawgn(EbNo_data_vec,'qam',max(4,M));
%  EsNo_vec=[EsNo_vec,EsNo_vec(end)+2,EsNo_vec(end)+3];
% 
% h_figure=figure(1);
% set(h_figure,'WindowStyle','Docked')
% 
% semilogy(EsNo_vec,BER_Ref,'-o')
% hold on
% grid on
% grid minor
% xlabel('EsNo [dB]')
% ylabel('BER')
% 
% 
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=LSChannel=Flat.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=LSChannel=MP 2 taps.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=MP 2 taps.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=LSChannel=MP 2 taps 18 chip 0.5.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=MP 2 taps 18 chip 0.5.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% 
% legend('Reference- Flat, No equalizer','Flat-Freq','Flat-Imp','MP-Freq',' MP-Imp','severe MP-Freq','severe MP-Imp')


% %%% MIMO Vs SISO
% 
% M=16;
% %load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=LSChannel=MP 2 taps 18 chip 0.5.mat')
%  load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat.mat')
% 
% 
% 
%  addition=2;
%  %EbNo_data_vec=[EbNo_data_vec,EbNo_data_vec(end)+2,EbNo_data_vec(end)+3];
%  BER_Ref = berawgn(EbNo_data_vec,'qam',max(4,M));
% % EsNo_vec=[EsNo_vec,EsNo_vec(end)+2,EsNo_vec(end)+3];
% 
% h_figure=figure(1);
% set(h_figure,'WindowStyle','Docked')
% 
% semilogy(EsNo_vec,BER_Ref,'-o')
% hold on
% grid on
% grid minor
% xlabel('EsNo [dB]')
% ylabel('BER')
% 
% 
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=2. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=3. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% 
% 
% 
% legend('Reference- SISO, No equalizer','SISO 1x1','MISO 1x2','MISO 1x3')

% %%% Fif: low versus high
% 
% M=16;
%  load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat. Fif=10.mat')
% 
% 
% 
%  addition=2;
%  %EbNo_data_vec=[EbNo_data_vec,EbNo_data_vec(end)+2,EbNo_data_vec(end)+3];
%  BER_Ref = berawgn(EbNo_data_vec,'qam',max(4,M));
% % EsNo_vec=[EsNo_vec,EsNo_vec(end)+2,EsNo_vec(end)+3];
% 
% h_figure=figure(1);
% set(h_figure,'WindowStyle','Docked')
% 
% semilogy(EsNo_vec,BER_Ref,'-o')
% hold on
% grid on
% grid minor
% xlabel('EsNo [dB]')
% ylabel('BER')
% 
% 
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat. Fif=10.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat. Fif=40.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat. Fif=40. inv sinc.mat')
% semilogy(EsNo_vec,BER_vec,'-o')
% hold on
% 
% 
% 
% legend('Reference- No equalizer, F_{IF}=10kHz','F_{IF}=10kHz','F_{IF}=40kHz','F_{IF}=40kHz. w inv sinc')
% 
% 


%%% Preambles enhancement

M=16;
%  load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=Flat. Fif=10.mat')
% 
% 
% 
%  addition=2;
%  %EbNo_data_vec=[EbNo_data_vec,EbNo_data_vec(end)+2,EbNo_data_vec(end)+3];
%  BER_Ref = berawgn(EbNo_data_vec,'qam',max(4,M));
% % EsNo_vec=[EsNo_vec,EsNo_vec(end)+2,EsNo_vec(end)+3];

h_figure=figure(1);
%set(h_figure,'WindowStyle','Docked')

% semilogy(EsNo_vec,BER_Ref,'-o')
% hold on
% grid on
% grid minor
% xlabel('EsNo [dB]')
% ylabel('BER')



load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=2. MIMO=1. Prmble enhance=0,0. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=MP.mat')
semilogy(EsNo_vec,BER_vec,'-o')
grid on
grid minor
xlabel('EsNo [dB]')
ylabel('BER')
hold on

load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=2. MIMO=1. Prmble enhance=3.65,2.8. Coding=0. Code rate=0.46667. Frequency offset=0.  Equalizer=MLChannel=MP.mat')
semilogy(EsNo_vec,BER_vec,'-o')
hold on



legend('Not Enhanced','Enhanced')




