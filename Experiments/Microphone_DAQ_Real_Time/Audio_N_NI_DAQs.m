%%%% Operarion of Audio devices to gether with NI's AFE devices

%%
close all
clear
clc

clear global BBB
global BBB

clear global CCC
global CCC

%% User Inputs

%%% Audio Devices

Fs_Audio=20e3*1; % minimum rate to sample audio well

N_frame_MIC=Fs_Audio*2;
N_frame_MIC=min(65536,2^(nextpow2(N_frame_MIC)));

%%% AFE

Fs_DAQ=5*Fs_Audio; % must no by pass AFE's capabilities, and on the other hand must be much higher than Audio devices to avoid gap (no real time action)


%%% NI configuration


ADC_Device='9215';
ADC_port=3;
ADC_FS=10;

DAC_port=0;
DAC_FS=10;


%% Audio objects definitions


Mic=dsp.AudioRecorder('SamplesPerfRame',N_frame_MIC,'SampleRate',Fs_Audio,...
    'NumChannels',1,'QueueDuration',1);


Speaker=dsp.AudioPlayer('SampleRate',Fs_Audio); % must have the same Fs as Microphone so that the voice does not distort

%% DAQ Sessions definitions

%%%%% Time period calculations: very important! (read comments)

Tspare=0.1;

T_Record_iteration_estimated_sec=N_frame_MIC*(1/Fs_Audio)
T_AFE_iteration_estimated_sec=N_frame_MIC*(1/Fs_DAQ) % based on emprical results; it takes the DAC always 0.2secs more then theory to empty its buffer

T_ADC_iteration_sec=T_AFE_iteration_estimated_sec+Tspare; % the ADC acquisition duration needs to be longer than DAC's, but shorter than the Audio device's
T_Comm_iteration_estimated_sec=T_ADC_iteration_sec+0.45 % empirical results show that AFE HW need 0.45 sec in addition to assigned comm time period to terminate their action completely and to enable recording from Mic again


if T_Comm_iteration_estimated_sec>T_Record_iteration_estimated_sec
    error('Fs_DAQ too low. Not enough room for transmission within 1 record frame')
end


%%% ADC 

ADC_Session = daq.createSession('ni');
ADC_Session.Rate = Fs_DAQ;
ADC_Session.DurationInSeconds =T_ADC_iteration_sec; % the ADC acquisition duration needs to be longer than DAC's, but shorter than the Audio device's



switch ADC_Device
    case '6212'
        
        Channel_ADC=ADC_Session.addAnalogInputChannel('Dev1',ADC_port,'Voltage');
    case '9215'
        
        ADC_port=strcat('ai',num2str(ADC_port));
        Channel_ADC=ADC_Session.addAnalogInputChannel('cDAQ1Mod1',ADC_port,'Voltage');
end

Channel_ADC.Range=[-ADC_FS,ADC_FS];
lh = addlistener(ADC_Session,'DataAvailable', @Store_Data1);

%%% DAC
DAC_Session = daq.createSession('ni');
DAC_Session.Rate = Fs_DAQ;

Channel_DAC=DAC_Session.addAnalogOutputChannel('Dev1', DAC_port, 'Voltage');


%% Recording and Communicating



%%%% Important: Mic object adapts itself automatically to the rest of the time delays
%%%% in the loop (AFE's time periods) so that total elapsed time is: N_frame_MIC*(1/Fs_Audio)
%%%% How? it records part of the time in the background (!) while other actions
%%%% are taking place (AFE communication mainly). As from 2nd iteration,
%%%% the AFE transmits on current iteration what had been recorded on
%%%% previous (!) iteration. E.g; say recording time is 2secs, and
%%%% required communication time including spare and AFE HW needs (0.45sec)
%%%% is 1.3 sec. Then, the microphone will record for 0.7sec frame #2,
%%%% from t=0 to t=1.3 with no other action taking place. at t=1.3 the AFE
%%%% will join in and transmit frame # 1 that had been recorded prviously.
%%%% from t=1.3 to t=2, both microphone anf AFE operate simultaneously but
%%%% on different data. At t=2, a new iteration begins; AFE shuts down and microphone begins to
%%%% record a new frame

n_max=5;
A=zeros(N_frame_MIC,1);
n=1;
t2=zeros(1,6);

disp('Sing into microphone now!')
disp('e_total should be equal to T_Record_iteration_estimated_sec as from 2nd iteration !')

while n<n_max
    
    %%%% Recording
    
    t1_record=clock;
    
    A=step(Mic); % acquisition of a frame of length N_frame_MIC
    
    t2_record=clock;
    e_record=etime(t2_record,t1_record); % time of only recording; follows is a simultaneous action of recording and AFE communicating
    
    %%% Communicating
    
    t1_comm=clock;
    
    A=A*9.5/max(abs(A)); %scaling to D/A's full scale
    queueOutputData(DAC_Session,A);
    BBB=0; % we need to zeroise that global variable so that it contains only current data acquisition form ADC
    startBackground(ADC_Session) % strat ADC acquisition
    startForeground(DAC_Session) % strat DAC transmission
    wait(ADC_Session)
    
%     %%%
%        t=0:1/Fs_Audio:(length(BBB)-1)*(1/Fs_Audio);
%         t=0:1/Fs_DAQ:(length(BBB)-1)*(1/Fs_DAQ);
%         plot(t/1e-3,BBB)
%         grid on;grid minor
%         xlabel('Time [msec]')
% %         t=0:1/Fs_Audio:(length(A)-1)*(1/Fs_Audio);
% %         hold on
% %         plot(t,A)
%     %%%
    
    t2_comm=clock;
    e_comm=etime(t2_comm,t1_comm);% the amount of time we miss because of non idelaity of AFE
    e_total=e_record+e_comm % should be equal more or less to T_Record_iteration_estimated_sec 

    if e_comm>e_record
        warning('Non Real Time') % we should have samples drops if that happens
    end
    
    n=n+1;
    
end

delete(lh)
release(Speaker)
release(Mic)


%% Display and send to speaker

DDD=CCC*(0.98/max(abs(CCC))); % speaker can only produce -1 to +1 volt

figure
t=0:1/Fs_Audio:(length(DDD)-1)*(1/Fs_Audio); % we want to display the acquired data in audio devices terms
set(gcf,'windowstyle','docked')
plot(t,DDD)
grid on;grid minor
title('Total Recorded Signal')
xlabel('Time [sec]')


step(Speaker,DDD) % if all goes well, the sequences will be a continuum of each other



