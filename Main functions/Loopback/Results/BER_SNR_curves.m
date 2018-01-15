
clc
clear
close all

M=16;
load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.95,2.95. Coding=1. Code rate=0.73333.mat')

addition=2;
EbNo_data_vec=[EbNo_data_vec,max(EbNo_data_vec)+addition];
BER_Ref = berawgn(EbNo_data_vec,'qam',max(4,M));
EsNo_vec=[EsNo_vec,EsNo_vec(end)+addition];


% BER_Ref = berawgn([EbNo_data_vec,max(EbNo_data_vec)+2],'qam',max(4,M));
% EsNo_vec=[EbNo_data_vec,max(EbNo_data_vec)+2]+EsNo_vec(1)-EbNo_data_vec(1);


h_figure=figure(1);
set(h_figure,'WindowStyle','Docked')
%subplot(1,2,1)
%h_axes_BER_EbNo=gca;
%h_line_BER_Ref_EbNo=semilogy([0:0.5:15],NaN(length([0:0.5:15]),1), 'g:o');
%hold on
%pause(1)
%h_line_BER_EbNo=semilogy([0:0.5:15],NaN(length([0:0.5:15]),1), ':o');
semilogy(EsNo_vec,BER_Ref,'-o')
hold on
grid on
grid minor
xlabel('EsNo [dB]')
ylabel('BER')
xlim([0,max(EsNo_vec)])
ylim(gca,[0.99e-5 1])
%legend('')

load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.95,2.95. Coding=1. Code rate=0.73333.mat')
semilogy(EsNo_vec,BER_vec,'-o')
hold on

load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.95,2.95. Coding=1. Code rate=0.6.mat')
semilogy(EsNo_vec,BER_vec,'-o')
hold on

load('N_FFT=256. N_CP=50. M_QAM=16. N Prmble long=8. MIMO=1. Prmble enhance=3.95,2.95. Coding=1. Code rate=0.46667.mat')
semilogy(EsNo_vec,BER_vec,'-o')
hold on


legend('no coding','coding rate=0.73','coding rate=0.6','coding rate=0.48')

% subplot(1,2,2)
% %set(h_figure,'WindowStyle','Docked')
% % h_axes_BER_EsNo=gca;
% % h_line_BER_Ref_EsNo=semilogy([0:0.5:15],NaN(length([0:0.5:15]),1), 'g:o');
% % hold on
% % pause(1)
% h_line_BER_EsNo=semilogy([0:0.5:15],NaN(length([0:0.5:15]),1), ':o');
% hold on
% grid on
% grid minor
% xlabel('EsNo [dB]')
% ylabel('BER')
% xlim([0,18])
% ylim(gca,[0.99e-5 1])