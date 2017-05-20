
close all
clear
clc


%%
dbstop if error
%dbstop if warning

clear persistent nn
clear global BBB

global BBB

%%


%Rx_PC
Tx_PC


%%
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
hCCDF = comm.CCDF('PAPROutputPort',true, 'MaximumPowerLimit', 2)  ;
[CCDFx,CCDFy,PAPR]=step(hCCDF,Signal_Tx_digital(Testing_data.Group_delay_Tx_total+1:end));
PAPR_dB=PAPR

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
plot(t_Rx/1e-3,Signal_Rx_digital)
xlabel('[msec]');grid on;grid minor
title(['Rx Signal at AFE interface.T preamble synch=',num2str(T_preamble_synch/1e-3),'[msec]. T preamble CE=',num2str(T_preamble_CE/1e-3),'[msec].'])

scatterplot(Symbol_stream_Rx);
grid on
grid minor
title(gca,['scatter plot. EVM=',num2str(EVM_dB),'[dB]'])
set(gcf,'windowstyle','docked');